/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "functions.h"


typedef struct
{
    Word16 inv_bin;
    Word16 numbytes;
    Word16 c_bp;
    Word16 c_bp_side;
    Word16 bytes;
    Word16 b_left;
    Word16 b_right;
    Word16 enc;
    Word16 bfi;
    Word16 be_bp_left;
    Word16 be_bp_right;
} Pc_State_fx;

typedef struct
{
    UWord32 ac_low_fx;
    UWord32 ac_range_fx;
    Word16  ac_cache_fx;
    Word16  ac_carry_fx;
    Word16  ac_carry_count_fx;
} Encoder_State_fx;

typedef struct
{
    UWord32 ac_low_fx;
    UWord32 ac_range_fx;
    UWord32 ac_help_fx;
    Word16  BER_detect;
    Pc_State_fx pc;
} Decoder_State_fx;

static void ac_dec_init_fx(UWord8 *ptr, Word16 *bp, Word16 *bp_side, Word16 *mask_side, Decoder_State_fx *st_fx /* i/o: Decoder State       */
);

static __forceinline void pc_init_fx(Word16 n_pc, Word16 numbytes, Word16 be_bp_left, Word16 be_bp_right, Word16 L_spec,
                                     Word16 enc, Word16 bfi, Pc_State_fx *pc /* i/o: Pc State */
);
static __forceinline Word16 check_pc_bytes(Word16 *bp, Word16 *bp_side, Word16 *mask_side, Word16 cur_bin, Word16 from_left,
                                           Pc_State_fx *pc /* i/o: Pc State */
);

static void ac_enc_init_fx(Encoder_State_fx *st_fx /* i/o: Encoder state       */
);

static void ac_enc_shift_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx /* i/o: Encoder state       */
);

static void write_indice_forward(UWord8 *ptr, Word16 bp, Word16 indice, Word16 numbits);

static void ac_encode_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx, /* i/o: Encoder state */
                         UWord32 cum_freq, /* i  : Cumulative frequency up to symbol   */
                         UWord32 sym_freq  /* i  : Symbol probability                  */
);

static Word16 ac_enc_finish_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx /* i/o: Encoder state       */
);

static Word16 ac_decode_fx(                         /* o  : Decoded cumulative frequency    */
                           Decoder_State_fx *st_fx, /* i/o: Decoder State                   */
                           Word16            pki);
static Word16 ac_decode_tns_order(                         /* o  : Decoded cumulative frequency    */
                                  Decoder_State_fx *st_fx, /* i/o: Decoder State                   */
                                  Word16            enable_lpc_weighting);
static Word16 ac_decode_tns_coef(                         /* o  : Decoded cumulative frequency    */
                                 Decoder_State_fx *st_fx, /* i/o: Decoder State                   */
                                 Word16            pki);
static Word16 ac_dec_update_fx(UWord8 *ptr, Word16 *bp, Word16 *bp_side, Word16 *mask_side, Word16 cur_bin,
                               Decoder_State_fx *st_fx, /* i/o: Decoder State           */
                               UWord32 cum_freq,        /* i  : Cumulative frequency    */
                               UWord32 sym_freq         /* i  : Symbol frequency        */
  );

/*************************************************************************/


Word16 processAriEncoder_fx(UWord8 *bytes, Word16 bp_side_in, Word16 mask_side_in, Word16 nbbits, Word16 xq[],
                            Word16 *tns_order, Word16 tns_numfilters, Word16 *tns_idx, Word16 lastnz,
                            Word16 *codingdata, Word8 *resBits, Word16 numResBits, Word16 lsbMode,
                            Word16 enable_lpc_weighting, Word8 *scratchBuffer)
{
    Dyn_Mem_Deluxe_In(
        Encoder_State_fx st;
        Word16           bp, bp_side, mask_side, extra_bits;
        Word16           a1, b1, a1_i, b1_i, a1_msb, b1_msb;
        Word16           lev1;
        Word16           nbits_side;
        Word16           tmp;
        Word16           fill_bits;
        UWord8 *         ptr;
        Word16           numResBitsEnc;
        Word16 *         lsb, nlsbs;
        Counter          i, n, k, lev;
    );

    lsb = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * MAX_LEN * 2 bytes */

    /* Init */
    a1_i          = 0; move16();
    b1_i          = 1; move16();
    bp            = 0; move16();
    numResBitsEnc = 0; move16();
    nlsbs         = 0; move16();
    ptr           = bytes;
    bp_side       = bp_side_in; move16();
    mask_side     = mask_side_in; move16();

    /*Start Encoding*/
    ac_enc_init_fx(&st);

    /* TNS data */
    FOR (n = 0; n < tns_numfilters; n++)
    {
        IF (tns_order[n] > 0)
        {
            ac_encode_fx(ptr, &bp, &st, ac_tns_order_cumfreq[enable_lpc_weighting][tns_order[n] - 1],
                         ac_tns_order_freq[enable_lpc_weighting][tns_order[n] - 1]);
            FOR (k = 0; k < tns_order[n]; k++)
            {
                ac_encode_fx(ptr, &bp, &st, ac_tns_coef_cumfreq[k][tns_idx[MAXLAG * n + k]],
                             ac_tns_coef_freq[k][tns_idx[MAXLAG * n + k]]);
            }
        }
    }

    IF (lsbMode == 0)
    {

        /*Main Loop through the 2-tuples*/
        FOR (k = 0; k < lastnz; k += 2)
        {
            IF (codingdata[1] < 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][0],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][0]);
            }
            ELSE IF (codingdata[1] == 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][codingdata[2]]);
                IF (xq[a1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[a1_i], 15));
                }
                IF (xq[b1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[b1_i], 15));
                }
            }
            ELSE IF (sub(codingdata[1], 1) == 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][VAL_ESC],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][VAL_ESC]);
                ac_encode_fx(ptr, &bp, &st,
                             ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]]);
                write_bit_backward(ptr, &bp_side, &mask_side, s_and(xq[a1_i], 1));
                write_bit_backward(ptr, &bp_side, &mask_side, s_and(xq[b1_i], 1));
                IF (xq[a1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[a1_i], 15));
                }
                IF (xq[b1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[b1_i], 15));
                }
            }
            ELSE
            {
                a1 = abs_s(xq[a1_i]);
                b1 = abs_s(xq[b1_i]);
                FOR (lev = 0; lev < codingdata[1]; lev++)
                {
                    lev1 = s_min(lev, 3);
                    ac_encode_fx(ptr, &bp, &st,
                                 ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC],
                                 ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC]);
                    write_bit_backward(ptr, &bp_side, &mask_side, s_and(shr_pos(a1, lev), 1));
                    write_bit_backward(ptr, &bp_side, &mask_side, s_and(shr_pos(b1, lev), 1));
                }
                lev1 = s_min(codingdata[1], 3);
                ac_encode_fx(ptr, &bp, &st,
                             ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]]);
                IF (xq[a1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[a1_i], 15));
                }
                IF (xq[b1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[b1_i], 15));
                }
            }

            a1_i += 2;
            b1_i += 2;
            codingdata += 3;

        } /*end of the 2-tuples loop*/
    }
    ELSE
    {
        /*Main Loop through the 2-tuples*/
        FOR (k = 0; k < lastnz; k += 2)
        {
            IF (codingdata[1] < 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][0],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][0]);
            }
            ELSE IF (codingdata[1] == 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][codingdata[2]]);
                IF (xq[a1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[a1_i], 15));
                }
                IF (xq[b1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[b1_i], 15));
                }
            }
            ELSE IF (sub(codingdata[1], 1) == 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][VAL_ESC],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][VAL_ESC]);
                ac_encode_fx(ptr, &bp, &st,
                             ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]]);
                a1_msb       = s_and(codingdata[2], 0x3);
                tmp          = s_and(xq[a1_i], 1);
                lsb[nlsbs++] = tmp; move16();
                test();
                IF (a1_msb == 0 && tmp > 0)
                {
                    if (xq[a1_i] > 0)
                    {
                        lsb[nlsbs++] = 0; move16();
                    }
                    if (xq[a1_i] < 0)
                    {
                        lsb[nlsbs++] = 1; move16();
                    }
                }
                IF (a1_msb != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[a1_i], 15));
                }
                b1_msb       = shr_pos(codingdata[2], 2);
                tmp          = s_and(xq[b1_i], 1);
                lsb[nlsbs++] = tmp; move16();
                test();
                IF (b1_msb == 0 && tmp > 0)
                {
                    if (xq[b1_i] > 0)
                    {
                        lsb[nlsbs++] = 0; move16();
                    }
                    if (xq[b1_i] < 0)
                    {
                        lsb[nlsbs++] = 1; move16();
                    }
                }
                IF (b1_msb != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[b1_i], 15));
                }
            }
            ELSE
            {
                a1           = abs_s(xq[a1_i]);
                b1           = abs_s(xq[b1_i]);
                a1_msb       = shr_pos(a1, 1);
                tmp          = s_and(a1, 1);
                lsb[nlsbs++] = tmp; move16();
                test();
                IF (a1_msb == 0 && tmp > 0)
                {
                    if (xq[a1_i] > 0)
                    {
                        lsb[nlsbs++] = 0; move16();
                    }
                    if (xq[a1_i] < 0)
                    {
                        lsb[nlsbs++] = 1; move16();
                    }
                }
                b1_msb       = shr_pos(b1, 1);
                tmp          = s_and(b1, 1);
                lsb[nlsbs++] = tmp; move16();
                test();
                IF (b1_msb == 0 && tmp > 0)
                {
                    if (xq[b1_i] > 0)
                    {
                        lsb[nlsbs++] = 0; move16();
                    }
                    if (xq[b1_i] < 0)
                    {
                        lsb[nlsbs++] = 1; move16();
                    }
                }
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[0]]][VAL_ESC],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[0]]][VAL_ESC]);
                FOR (lev = 1; lev < codingdata[1]; lev++)
                {
                    lev1 = s_min(lev, 3);
                    ac_encode_fx(ptr, &bp, &st,
                                 ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC],
                                 ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC]);
                    write_bit_backward(ptr, &bp_side, &mask_side, s_and(shr_pos(a1, lev), 1));
                    write_bit_backward(ptr, &bp_side, &mask_side, s_and(shr_pos(b1, lev), 1));
                }
                lev1 = s_min(codingdata[1], 3);
                ac_encode_fx(ptr, &bp, &st,
                             ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]]);
                IF (a1_msb != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[a1_i], 15));
                }
                IF (b1_msb != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[b1_i], 15));
                }
            }

            a1_i += 2;
            b1_i += 2;
            codingdata += 3;

        } /*end of the 2-tuples loop*/
    }

    /* Side bits (in sync with the decoder) */
    nbits_side = sub(nbbits, add(shl_pos(bp_side, 3), sub(norm_s(mask_side), 6)));

    /* Residual bits (in sync with the decoder) */
    extra_bits = sub(norm_ul(st.ac_range_fx), 6);
    if (st.ac_cache_fx >= 0)
    {
        extra_bits = add(extra_bits, 8);
    }
    if (st.ac_carry_count_fx > 0)
    {
        extra_bits = add(extra_bits, shl_pos(st.ac_carry_count_fx, 3));
    }
    n = sub(nbbits, add(shl_pos(bp, 3), add(extra_bits, nbits_side))); move16();
    assert(n >= 0);

    IF (lsbMode == 0)
    {
        numResBitsEnc = s_min(numResBits, n);
        FOR (i = 0; i < numResBitsEnc; i++)
        {
            write_bit_backward(ptr, &bp_side, &mask_side, (Word16)resBits[i]);
        }
    }
    ELSE
    {
        nlsbs = s_min(nlsbs, n);
        FOR (k = 0; k < nlsbs; k++)
        {
            write_bit_backward(ptr, &bp_side, &mask_side, lsb[k]);
        }
    }

    /* End arithmetic coder, overflow management */
    extra_bits = ac_enc_finish_fx(ptr, &bp, &st);

    /* Fill bits (for debugging, the exact number of fill bits cannot be computed in the decoder)*/
    fill_bits = nbbits - (bp * 8 + extra_bits + nbits_side + nlsbs + numResBitsEnc);

    Dyn_Mem_Deluxe_Out();
    return fill_bits;
}


void processAriDecoder_fx(UWord8 *bytes, Word16 *bp_side, Word16 *mask_side, Word16 nbbits, Word16 L_spec,
                          Word16 fs_idx, Word16 enable_lpc_weighting, Word16 tns_numfilters, Word16 lsbMode,
                          Word16 lastnz, Word16 *bfi, Word16 *tns_order, Word16 fac_ns_idx, Word16 gg_idx,
                          Word16 frame_dms,
                          Word16 n_pc, Word16 be_bp_left, Word16 be_bp_right, Word16 enc, Word16 *spec_inv_idx, Word16 *b_left,
                          Word16 *resBits, Word16 *x, Word16 *nf_seed, Word16 *resQdata, Word16 *tns_idx,
                          Word16 *zero_frame, Word8 *scratchBuffer)
{

    Decoder_State_fx st;
    Word16  a, b, t, a1, b1, a1_i, b1_i, bp;
    Word16  esc_nb;
    Word16  rateFlag;
    Word16  r;
    Word16  nt_half;
    Word16  c;
    Word16  nbits_side, extra_bits, nbits_ari;
    UWord8 *ptr;
    Word32  tmp32;
    Word16  lsb_ind_c;
    Word16 *lsb_ind;
    Word16  tmp;
    Counter n, k, lev;
    Counter i;

#ifdef DYNMEM_COUNT
    struct _dynmem
    {
        Decoder_State_fx st;
        Pc_State_fx pc;

        Word16  a, b, t, a1, b1, a1_i, b1_i, bp;
        Word16  esc_nb;
        Word16  rateFlag;
        Word16  r;
        Word16  nt_half;
        Word16  c;
        Word16  nbits_side, extra_bits, nbits_ari;
        UWord8 *ptr;
        Word32  tmp32;
        Word16  lsb_ind_c;
        Word16 *lsb_ind;
        Word16  tmp;
        Counter i, n, k, lev;
    };
    Dyn_Mem_In("processAriDecoder_fx", sizeof(struct _dynmem));
#endif


    lsb_ind = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size 2 * MAX_LEN bytes */

    /* Rate flag */
    rateFlag = 0; move16();
    if (sub(nbbits, add(160, i_mult(fs_idx, 160))) > 0)
    {
        rateFlag = 2 << NBITS_CONTEXT; move16();
    }

    pc_init_fx(n_pc, shr_pos(nbbits,3), be_bp_left, be_bp_right, L_spec, enc, *bfi, &st.pc);

    /* Init */
    nt_half   = shr_pos(L_spec, 1);
    c         = 0; move16();
    t         = 0; move16();
    a1_i      = 0; move16();
    b1_i      = 1; move16();
    bp        = 0; move16();
    if (enc == 0)
    {
        bp = add(bp, st.pc.bytes);  move16();
    }
    *spec_inv_idx = L_spec; move16();
    *b_left = -1; move16();
    lsb_ind_c = 0; move16();

    ptr = bytes;

    /* Start Decoding */
    ac_dec_init_fx(ptr, &bp, bp_side, mask_side, &st);

    /* Decode TNS data */
#ifdef NONBE_BER_DETECT
    tmp = MAXLAG;
    IF (sub(frame_dms, 25) == 0)
    {
        tmp = shr_pos(tmp, 1);
    }
    IF (sub(frame_dms, 50) == 0)
    {
        tmp = shr_pos(tmp, 1);
    }
#endif
    FOR (n = 0; n < tns_numfilters; n++)
    {
        IF (tns_order[n] > 0)
        {
            tns_order[n] = ac_decode_tns_order(&st, enable_lpc_weighting); move16();
            tns_order[n] = add(tns_order[n], 1);                           move16();
#ifdef NONBE_BER_DETECT
            IF (tns_order[n] > tmp)
            {
                GOTO ber_detect;
            }
#endif
            if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, 0, &st, ac_tns_order_cumfreq[enable_lpc_weighting][tns_order[n] - 1],
                             ac_tns_order_freq[enable_lpc_weighting][tns_order[n] - 1]) != 0)
            {
                GOTO ber_detect;
            }
            FOR (k = 0; k < tns_order[n]; k++)
            {
                IF (sub(*bp_side, bp) < 0)
                {
                     GOTO ber_detect;
                }
                tns_idx[MAXLAG * n + k] = ac_decode_tns_coef(&st, k); move16();
                if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, 0, &st, ac_tns_coef_cumfreq[k][tns_idx[MAXLAG * n + k]],
                                 ac_tns_coef_freq[k][tns_idx[MAXLAG * n + k]]) != 0)
                {
                    GOTO ber_detect;
                }
            }
        }
    }
    IF (st.BER_detect > 0)
    {
        GOTO ber_detect;
    }

    IF (lsbMode == 0)
    {

        /*Main Loop through the 2-tuples*/
        FOR (k = 0; k < lastnz; k += 2)
        {

            /* Get context */
            t = add(c, rateFlag);
            if (sub(k, nt_half) > 0)
            {
                t = add(t, 1 << NBITS_CONTEXT);
            }

            r = ac_decode_fx(&st, ari_spec_lookup[t]);
            if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st, ari_spec_cumfreq[ari_spec_lookup[t]][r],
                             ari_spec_freq[ari_spec_lookup[t]][r]) != 0)
            {
                GOTO ber_detect;
            }

            IF (r == 0)
            {
                x[a1_i] = 0; move16();
                x[b1_i] = 0; move16();
                c       = add(shl_pos(s_and(c, 0xf), 4), 1);
            }
            ELSE IF (sub(r, VAL_ESC) < 0)
            {
                a = s_and(r, 0x3);
                b = shr_pos(r, 2);
                c = add(shl_pos(s_and(c, 0xf), 4), add(add(a, b), 1));
                IF (a > 0)
                {
                    if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                    {
                        GOTO ber_detect;
                    }
                    if (read_bit(ptr, bp_side, mask_side) != 0)
                    {
                        a = negate(a);
                    }
                }
                x[a1_i] = a; move16();
                IF (b > 0)
                {
                    if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0)
                    {
                        GOTO ber_detect;
                    }
                    if (read_bit(ptr, bp_side, mask_side) != 0)
                    {
                        b = negate(b);
                    }
                }
                x[b1_i] = b; move16();
            }
            ELSE
            {
                if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                {
                    GOTO ber_detect;
                }
                a = read_bit(ptr, bp_side, mask_side);
                if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                {
                    GOTO ber_detect;
                }
                b = read_bit(ptr, bp_side, mask_side);
                r = ac_decode_fx(&st, ari_spec_lookup[t + Tab_esc_nb[1]]);
                if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st, ari_spec_cumfreq[ari_spec_lookup[t + Tab_esc_nb[1]]][r],
                                 ari_spec_freq[ari_spec_lookup[t + Tab_esc_nb[1]]][r]) != 0)
                {
                    GOTO ber_detect;
                }
                IF (sub(r, VAL_ESC) < 0)
                {
                    a1 = s_and(r, 0x3);
                    b1 = shr_pos(r, 2);
                    a  = add(shl_pos(a1, 1), a);
                    b  = add(shl_pos(b1, 1), b);
                    IF (a > 0)
                    {
                        if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
                            a = negate(a);
                        }
                    }
                    x[a1_i] = a; move16();
                    IF (b > 0)
                    {
                        if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
                            b = negate(b);
                        }
                    }
                    x[b1_i] = b; move16();
                    c       = add(shl_pos(s_and(c, 0xf), 4), add(shl_pos(add(a1, b1), 1), 1));
                }
                ELSE
                {
                    if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                    {
                        GOTO ber_detect;
                    }
                    a = add(shl_pos(read_bit(ptr, bp_side, mask_side), 1), a);
                    if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                    {
                        GOTO ber_detect;
                    }
                    b = add(shl_pos(read_bit(ptr, bp_side, mask_side), 1), b);
                    FOR (lev = 2; lev < 14; lev++)
                    {
                        esc_nb = s_min(lev, 3);
                        r      = ac_decode_fx(&st, ari_spec_lookup[t + Tab_esc_nb[esc_nb]]);
                        if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st, ari_spec_cumfreq[ari_spec_lookup[t + Tab_esc_nb[esc_nb]]][r],
                                         ari_spec_freq[ari_spec_lookup[t + Tab_esc_nb[esc_nb]]][r]) != 0)
                        {
                            GOTO ber_detect;
                        }
                        IF (sub(r, VAL_ESC) < 0)
                        {
                            BREAK;
                        }
                        if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }
                        a = add(shl(read_bit(ptr, bp_side, mask_side), lev), a);
                        if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }
                        b = add(shl(read_bit(ptr, bp_side, mask_side), lev), b);
                    }
                    /* check for bitflip */
                    IF (sub(lev, 14) == 0)
                    {
                        GOTO ber_detect;
                    }

                    b1 = shr_pos(r, 2);
                    a1 = s_and(r, 0x3);
                    a  = add(shl(a1, lev), a);
                    b  = add(shl(b1, lev), b);
                    IF (a > 0)
                    {
                        if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
                            a = negate(a);
                        }
                    }
                    x[a1_i] = a; move16();
                    IF (b > 0)
                    {
                        if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
                            b = negate(b);
                        }
                    }
                    x[b1_i] = b; move16();
                    c       = add(shl_pos(s_and(c, 0xf), 4), add(esc_nb, 12));
                }
            }

            test();test();
            IF ((sub(sub(bp, *bp_side), 3) > 0 && sub(st.pc.c_bp, st.pc.c_bp_side) == 0) || st.BER_detect > 0)
            {
                GOTO ber_detect;
            }

            a1_i += 2;
            b1_i += 2;
        }
    }
    ELSE
    {
        /*Main Loop through the 2-tuples*/
        FOR (k = 0; k < lastnz; k += 2)
        {

            /* Get context */
            t = add(c, rateFlag);
            if (sub(k, nt_half) > 0)
            {
                t = add(t, 1 << NBITS_CONTEXT);
            }

            r = ac_decode_fx(&st, ari_spec_lookup[t]);
            if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st, ari_spec_cumfreq[ari_spec_lookup[t]][r],
                             ari_spec_freq[ari_spec_lookup[t]][r]) != 0)
            {
                GOTO ber_detect;
            }

            IF (r == 0)
            {
                x[a1_i] = 0; move16();
                x[b1_i] = 0; move16();
                c       = add(shl_pos(s_and(c, 0xf), 4), 1);
            }
            ELSE IF (sub(r, VAL_ESC) < 0)
            {
                a = s_and(r, 0x3);
                b = shr_pos(r, 2);
                c = add(shl_pos(s_and(c, 0xf), 4), add(add(a, b), 1));
                IF (a > 0)
                {
                    if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                    {
                        GOTO ber_detect;
                    }
                    if (read_bit(ptr, bp_side, mask_side) != 0)
                    {
                        a = negate(a);
                    }
                }
                x[a1_i] = a; move16();
                IF (b > 0)
                {
                    if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0)
                    {
                        GOTO ber_detect;
                    }
                    if (read_bit(ptr, bp_side, mask_side) != 0)
                    {
                        b = negate(b);
                    }
                }
                x[b1_i] = b; move16();
            }
            ELSE
            {
                r = ac_decode_fx(&st, ari_spec_lookup[t + Tab_esc_nb[1]]);
                if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st, ari_spec_cumfreq[ari_spec_lookup[t + Tab_esc_nb[1]]][r],
                                 ari_spec_freq[ari_spec_lookup[t + Tab_esc_nb[1]]][r]) != 0)
                {
                    GOTO ber_detect;
                }
                IF (sub(r, VAL_ESC) < 0)
                {
                    a1 = s_and(r, 0x3);
                    b1 = shr_pos(r, 2);
                    a  = shl_pos(a1, 1);
                    b  = shl_pos(b1, 1);
                    IF (a > 0)
                    {
                        if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
                            a = negate(a);
                        }
                    }
                    x[a1_i] = a; move16();
                    IF (b > 0)
                    {
                        if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
                            b = negate(b);
                        }
                    }
                    x[b1_i]              = b; move16();
                    c                    = add(shl_pos(s_and(c, 0xf), 4), add(shl_pos(add(a1, b1), 1), 1));
                    lsb_ind[lsb_ind_c++] = k; move16();
                }
                ELSE
                {
                    if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                    {
                        GOTO ber_detect;
                    }
                    a = shl_pos(read_bit(ptr, bp_side, mask_side), 1);
                    if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                    {
                        GOTO ber_detect;
                    }
                    b = shl_pos(read_bit(ptr, bp_side, mask_side), 1);
                    FOR (lev = 2; lev < 14; lev++)
                    {
                        esc_nb = s_min(lev, 3);
                        r      = ac_decode_fx(&st, ari_spec_lookup[t + Tab_esc_nb[esc_nb]]);
                        if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st, ari_spec_cumfreq[ari_spec_lookup[t + Tab_esc_nb[esc_nb]]][r],
                                         ari_spec_freq[ari_spec_lookup[t + Tab_esc_nb[esc_nb]]][r]) != 0)
                        {
                            GOTO ber_detect;
                        }
                        IF (sub(r, VAL_ESC) < 0)
                        {
                            BREAK;
                        }
                        if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }
                        a = add(shl(read_bit(ptr, bp_side, mask_side), lev), a);
                        if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }
                        b = add(shl(read_bit(ptr, bp_side, mask_side), lev), b);
                    }
                    /* check for bitflip */
                    IF (sub(lev, 14) == 0)
                    {
                        GOTO ber_detect;
                    }

                    b1 = shr_pos(r, 2);
                    a1 = s_and(r, 0x3);
                    a  = add(shl(a1, lev), a);
                    b  = add(shl(b1, lev), b);
                    IF (a > 0)
                    {
                        if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
                            a = negate(a);
                        }
                    }
                    x[a1_i] = a; move16();
                    IF (b > 0)
                    {
                        if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
                            b = negate(b);
                        }
                    }
                    x[b1_i]              = b; move16();
                    c                    = add(shl_pos(s_and(c, 0xf), 4), add(esc_nb, 12));
                    lsb_ind[lsb_ind_c++] = k; move16();
                }
            }

            test();test();
            IF ((sub(sub(bp, *bp_side), 3) > 0 && sub(st.pc.c_bp, st.pc.c_bp_side) == 0) || st.BER_detect > 0)
            {
                GOTO ber_detect;
            }

            a1_i += 2;
            b1_i += 2;
        }
    }

    IF (L_spec > k)
    {
        basop_memset(&x[k], 0, (L_spec - k) * sizeof(*x));
    }

    nbits_side = sub(nbbits, add(shl_pos(*bp_side, 3), sub(norm_s(*mask_side), 6)));
    extra_bits  = sub(norm_ul(st.ac_range_fx), 6);
    nbits_ari   = shl_pos(sub(bp, 3), 3);
    IF (enc == 0)
    {
        IF (st.pc.c_bp == 0)
        {
            nbits_ari = shl_pos(sub(sub(bp, st.pc.bytes), 3), 3);
        }
        ELSE
        {
            nbits_ari = shl_pos(add(bp, sub(sub(st.pc.b_left, st.pc.bytes), 3)), 3);
        }

        IF (st.pc.c_bp_side != 0)
        {
            nbits_side = sub(add(sub(nbbits, shl_pos(st.pc.b_left, 3)), shl_pos(sub(st.pc.bytes, *bp_side), 3)), sub(norm_s(*mask_side), 6));
        }
    }

    n = sub(nbbits, add(nbits_ari, add(extra_bits, nbits_side))); move16();

    IF (n < 0)
    {
        GOTO ber_detect;
    }

    IF (lsbMode == 0)
    {
        *resBits = n; move16();
        FOR (k = 0; k < L_spec; k++)
        {
            IF (x[k] != 0)
            {
                IF (n == 0)
                {
                    BREAK;
                }
                if (check_pc_bytes(&bp, bp_side, mask_side, st.pc.inv_bin, 0, &st.pc) != 0)
                {
                    GOTO ber_detect_res;
                }
                *resQdata++ = read_bit(ptr, bp_side, mask_side); move16();
                n           = sub(n, 1);
            }
        }
        *resBits = sub(*resBits, n);
    }
    ELSE
    {
        *resBits = 0;
        FOR (k = 0; k < lsb_ind_c; k++)
        {
            a = x[lsb_ind[k]]; move16();
            IF (n == 0)
            {
                BREAK;
            }
            if (check_pc_bytes(&bp, bp_side, mask_side, st.pc.inv_bin, 0, &st.pc) != 0)
            {
                GOTO ber_detect_res;
            }
            tmp = read_bit(ptr, bp_side, mask_side);
            n   = sub(n, 1);
            IF (tmp > 0)
            {
                if (a > 0)
                {
                    a = add(a, 1);
                }
                if (a < 0)
                {
                    a = sub(a, 1);
                }
                IF (a == 0)
                {
                    IF (n == 0)
                    {
                        BREAK;
                    }
                    a = 1;
                    if (check_pc_bytes(&bp, bp_side, mask_side, st.pc.inv_bin, 0, &st.pc) != 0)
                    {
                        GOTO ber_detect_res;
                    }
                    if (read_bit(ptr, bp_side, mask_side) != 0)
                    {
                        a = negate(a);
                    }
                    n = sub(n, 1);
                }
            }

            x[lsb_ind[k]] = a;                 move16();
            b             = x[lsb_ind[k] + 1]; move16();
            IF (n == 0)
            {
                BREAK;
            }
            if (check_pc_bytes(&bp, bp_side, mask_side, st.pc.inv_bin, 0, &st.pc) != 0)
            {
                GOTO ber_detect_res;
            }
            tmp = read_bit(ptr, bp_side, mask_side);
            n   = sub(n, 1);
            IF (tmp > 0)
            {
                if (b > 0)
                {
                    b = add(b, 1);
                }
                if (b < 0)
                {
                    b = sub(b, 1);
                }
                IF (b == 0)
                {
                    IF (n == 0)
                    {
                        BREAK;
                    }
                    b = 1;
                    if (check_pc_bytes(&bp, bp_side, mask_side, st.pc.inv_bin, 0, &st.pc) != 0)
                    {
                        GOTO ber_detect_res;
                    }
                    if (read_bit(ptr, bp_side, mask_side) != 0)
                    {
                        b = negate(b);
                    }
                    n = sub(n, 1);
                }
            }
            x[lsb_ind[k] + 1] = b; move16();
        }
    }

    /* Noise Filling seed */
    tmp32 = L_deposit_l(0);
    FOR (i = 0; i < L_spec; i++)
    {
        tmp32 = L_mac0(tmp32, abs_s(x[i]), i);
    }
    *nf_seed = extract_l(tmp32); move16();

    /* Detect zero frame */
    test(); test(); test(); test();
    IF (sub(lastnz, 2) == 0 && sub(x[0], 0) == 0 && sub(x[1], 0) == 0 && sub(gg_idx, 0) == 0 && sub(fac_ns_idx, 7) == 0)
    {
        *zero_frame = 1; move16();
    }
    ELSE
    {
        *zero_frame = 0; move16();
    }

    IF (enc)
    {
        IF (st.pc.bytes > 0)
        {
            IF (sub(st.pc.b_left, shr_pos(nbbits,3)) > 0)
            {
                *b_left = sub(*bp_side, st.pc.bytes);
            }
        }
    }

    IF (sub(*bfi, 2) == 0)
    {
        IF (sub(*spec_inv_idx, L_spec) == 0)
        {
            *bfi = 0;
        }
    }
    GOTO bail;

/* goto for bit error handling */
ber_detect:
    *bfi = 1; move16();
    *b_left = st.pc.b_left; move16();
    test();
    IF (st.pc.inv_bin > 0 && sub(st.pc.inv_bin, L_spec) <= 0)
    {
        *spec_inv_idx = st.pc.inv_bin; move16();
        *bfi          = 2; move16();
        *resBits      = 0; move16();
        *zero_frame   = 0; move16();
        /* Noise Filling seed */
        tmp32 = L_deposit_l(0);
        FOR (i = 0; i < *spec_inv_idx; i++)
        {
            tmp32 = L_mac0(tmp32, abs_s(x[i]), i);
        }
        *nf_seed = extract_l(tmp32); move16();
    }
    GOTO bail;

/* goto for bit error handling in residual signal */
ber_detect_res:
    *b_left     = st.pc.b_left; move16();
    *resBits    = 0; move16();
    *bfi        = 0; move16();
    *zero_frame = 0; move16();
    /* Noise Filling seed */
    tmp32 = L_deposit_l(0);
    FOR (i = 0; i < *spec_inv_idx; i++)
    {
        tmp32 = L_mac0(tmp32, abs_s(x[i]), i);
    }
    *nf_seed = extract_l(tmp32); move16();
    GOTO bail;
    
    /* goto, because of dynmem out */
bail:
    Dyn_Mem_Deluxe_Out();
}


void processAriDecoderScaling_fx(Word16 *data16, Word16 dataLen, Word32 *data32, Word16 *data_e)
{
    Counter i;
    Dyn_Mem_Deluxe_In(
        Word16 tmp, shift;
        Word16 x_min, x_max;
    );

    x_max = 0; move16();
    x_min = 0; move16();

    FOR (i = 0; i < dataLen; i++)
    {
        if (data16[i] > 0)
            x_max = s_max(x_max, data16[i]);
        if (data16[i] < 0)
            x_min = s_min(x_min, data16[i]);
    }

    tmp   = s_max(x_max, negate(x_min));
    shift = norm_s(tmp);
    if (tmp == 0)
    {
        shift = 15; move16();
    }

    FOR (i = 0; i < dataLen; i++)
    {
        data32[i] = L_shl_pos(L_deposit_h(data16[i]), shift); move32();
    }

    *data_e = sub(15, shift); move16();

    Dyn_Mem_Deluxe_Out();
}


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

/*************************************************************************/

static __forceinline UWord32 UL_addNs24(UWord32 UL_var1, UWord32 UL_var2, UWord16 *wrap)
{
    return UL_lshr(UL_addNs(UL_lshl(UL_var1, 8), UL_lshl(UL_var2, 8), wrap), 8);
}

Word16 find_last_nz_pair(const Word16 x[], Word16 length)
{
    Dyn_Mem_Deluxe_In(
        Word16  last_nz, lobs[4];
        Counter stage, i;
    );

    lobs[0] = 4;                  move16();
    lobs[1] = shr_pos(length, 1); /* length/2 */
    move16();
    lobs[2] = add(lobs[1], shr_pos(length, 2)); move16();
    lobs[3] = add(lobs[2], shr_pos(length, 3)); move16();

    last_nz = 0;      move16();
    i       = length; move16();
    FOR (stage = 3; stage >= 0; --stage)
    {
        /* unmapped kernel */
        FOR (; i >= lobs[stage]; i -= 2)
        {
            if (x[i - 2] != 0)
            {
                last_nz = s_max(last_nz, i);
            }
            if (x[i - 1] != 0)
            {
                last_nz = s_max(last_nz, i);
            }
        }
        IF (last_nz > 0)
        {
            BREAK;
        }
    }

    Dyn_Mem_Deluxe_Out();
    return s_max(last_nz, 2);
}


void write_bit_backward(UWord8 *ptr, Word16 *bp, Word16 *mask, Word16 bit)
{
    if (bit > 0)
    {
        ptr[*bp] = (UWord8)s_or((Word16)ptr[*bp], *mask); move16();
    }
    *mask = lshl_pos(*mask, 1); move16();
    if (sub(*mask, 0x100) == 0)
    {
        *mask = 1; move16();
    }
    if (sub(*mask, 1) == 0)
    {
        *bp = sub(*bp, 1); move16();
    }
}


void write_indice_backward(UWord8 *ptr, Word16 *bp, Word16 *mask, Word16 indice, Word16 numbits)
{
    Dyn_Mem_Deluxe_In(
        Counter k;
        Word16 bit;
    );

    FOR (k = 0; k < numbits; k++)
    {
        bit = s_and(indice, 1);
        write_bit_backward(ptr, bp, mask, bit);
        indice = lshr(indice, 1);
    }

    Dyn_Mem_Deluxe_Out();
}


static __forceinline void write_indice_forward(UWord8 *ptr, Word16 bp, Word16 indice, Word16 numbits)
{
    Dyn_Mem_Deluxe_In(
        Counter k;
        Word16  bit, mask, tmp;
    );

    tmp  = (Word16)ptr[bp]; move16();
    mask = 0x80;            move16();
    FOR (k = 0; k < numbits; k++)
    {
        bit = s_and(indice, mask);
        tmp = s_or(tmp, mask);
        if (bit == 0)
        {
            tmp = sub(tmp, mask);
        }
        mask = lshr(mask, 1);
    }
    ptr[bp] = (UWord8)tmp; move16();

    Dyn_Mem_Deluxe_Out();
}

static __forceinline void ac_enc_init_fx(Encoder_State_fx *st_fx) /* i/o: Encoder state       */
{
    st_fx->ac_low_fx         = L_deposit_l(0); move32();
    st_fx->ac_range_fx       = 0x00ffffff;     move32();
    st_fx->ac_cache_fx       = -1;             move16();
    st_fx->ac_carry_fx       = 0;              move16();
    st_fx->ac_carry_count_fx = 0;              move16();
}

static __forceinline void ac_enc_shift_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx) /* i/o: Encoder state */
{
    test();
    L_sub(0, 0); /* For comparision in if */
    IF (st_fx->ac_low_fx < (0x00ff0000UL) || sub(st_fx->ac_carry_fx, 1) == 0)
    {
        IF (st_fx->ac_cache_fx >= 0)
        {
            ptr[(*bp)++] = (UWord8)add(st_fx->ac_cache_fx, st_fx->ac_carry_fx); move16();
        }

        WHILE (st_fx->ac_carry_count_fx > 0)
        {
            ptr[(*bp)++]             = (UWord8)s_and(add(st_fx->ac_carry_fx, 0xff), 255); move16();
            st_fx->ac_carry_count_fx = sub(st_fx->ac_carry_count_fx, 1);                  move16();
        }

        st_fx->ac_cache_fx = u_extract_l(UL_lshr_pos(st_fx->ac_low_fx, 16)); move16();
        st_fx->ac_carry_fx = 0;                                              move16();
    }
    ELSE
    {
        st_fx->ac_carry_count_fx = add(st_fx->ac_carry_count_fx, 1); move16();
    }
    st_fx->ac_low_fx = UL_and(UL_lshl_pos(st_fx->ac_low_fx, 8), 0x00ffffff); move32();
}

static __forceinline void ac_encode_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx, /* i/o: Encoder state */
                                       UWord32 cum_freq, /* i  : Cumulative frequency up to symbol   */
                                       UWord32 sym_freq)  /* i  : Symbol probability                  */
{
    Dyn_Mem_Deluxe_In(
        UWord32 r, tmp;
        UWord16 carry;
    );

    r   = UL_lshr_pos(st_fx->ac_range_fx, 10);
    tmp = UL_Mpy_32_32(r, cum_freq);

    assert(r < (1U << 24));
    assert(cum_freq < (1U << 24));
    assert(tmp < (1U << 24));
    assert(st_fx->ac_low_fx < (1U << 24));
    st_fx->ac_low_fx = UL_addNs24(st_fx->ac_low_fx, tmp, &carry); move32();

    if (carry != 0)
    {
        st_fx->ac_carry_fx = carry; move16();
    }

    st_fx->ac_range_fx = UL_Mpy_32_32(r, sym_freq); move32();

    assert(cum_freq < (1U << 24));
    assert(st_fx->ac_range_fx < (1U << 24));
    WHILE (st_fx->ac_range_fx < (1U << 16))
    {
        L_sub(0, 0); /* Comparison in while */
        st_fx->ac_range_fx = UL_lshl_pos(st_fx->ac_range_fx, 8); move32();

        assert(st_fx->ac_range_fx < (1U << 24));

        ac_enc_shift_fx(ptr, bp, st_fx);
    }

    Dyn_Mem_Deluxe_Out();
}

static __forceinline Word16 ac_enc_finish_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx) /* i/o: Encoder state */
{
    Dyn_Mem_Deluxe_In(
        UWord32 val, mask, high;
        Word16  bits;
        UWord16 over1, over2;
    );

    /*bits = 24 - log2_i(st->ac_range); */
    bits = sub(norm_ul(st_fx->ac_range_fx), 7);

    mask = UL_lshr(0x00ffffff, bits);

    val  = UL_addNs24(st_fx->ac_low_fx, mask, &over1);
    high = UL_addNs24(st_fx->ac_low_fx, st_fx->ac_range_fx, &over2);

    L_xor(0, 0);    /* For bit not */
    UL_and(1U, 1U); /* added counters */
    val = L_and(val, (~mask) & 0x00ffffff);

    L_xor(0, 0); /* For bit not */
    IF ((L_xor(over1, over2)) == 0)
    {
        L_sub(0, 0); /* For comparision in if */
        IF (UL_addNsD(val, mask) >= high)
        {
            bits = add(bits, 1);
            mask = UL_lshr_pos(mask, 1);
            val  = UL_and(UL_addNsD(st_fx->ac_low_fx, mask), (~mask) & 0x00ffffff);
            L_xor(0, 0);
            UL_and(1, 1); /* For bit not , mask */
        }

        if (val < st_fx->ac_low_fx)
        {
            st_fx->ac_carry_fx = 1; move16();
        }
    }

    st_fx->ac_low_fx = val; move32();

    FOR (; bits > 0; bits -= 8)
    {
        ac_enc_shift_fx(ptr, bp, st_fx);
    }
    bits = add(bits, 8);

    assert(st_fx->ac_carry_fx == 0);

    IF (st_fx->ac_carry_count_fx > 0)
    {
        ptr[(*bp)++] = (UWord8)st_fx->ac_cache_fx; move16();

        FOR (; st_fx->ac_carry_count_fx > 1; st_fx->ac_carry_count_fx--)
        {
            ptr[(*bp)++] = 0xff; move16();
        }
        write_indice_forward(ptr, *bp, lshr(0xff, sub(8, bits)), bits);
    }
    ELSE
    {
        write_indice_forward(ptr, *bp, st_fx->ac_cache_fx, bits);
    }

    Dyn_Mem_Deluxe_Out();
    return bits;
}


__forceinline Word16 read_bit(UWord8 *ptr, Word16 *bp, Word16 *mask)
{
    Dyn_Mem_Deluxe_In(
        Word16 bit;
    );

    bit = 0; move16();
    if (s_and((Word16)ptr[*bp], *mask) > 0)
    {
        bit = 1; move16();
    }
    *mask = lshl_pos(*mask, 1); move16();
    if (sub(*mask, 0x100) == 0)
    {
        *mask = 1; move16();
    }
    if (sub(*mask, 1) == 0)
    {
        *bp = sub(*bp, 1); move16();
    }


    Dyn_Mem_Deluxe_Out();
    return bit;
}


static __forceinline void ac_dec_init_fx(UWord8 *ptr, Word16 *bp, Word16 *bp_side, Word16 *mask_side, Decoder_State_fx *st_fx) /* i/o: Decoder State */
{
    Dyn_Mem_Deluxe_In(
        Counter i;
    );


    st_fx->ac_low_fx = L_deposit_l(0); move32();

    st_fx->ac_range_fx = 0x00ffffff; move32();
    FOR (i = 0; i < 3; i++)
    {
        if (check_pc_bytes(bp, bp_side, mask_side, 0, 1, &st_fx->pc) != 0)
        {
            Dyn_Mem_Deluxe_Out();
            return;
        }
        st_fx->ac_low_fx = UL_addNsD(UL_lshl_pos(st_fx->ac_low_fx, 8), UL_deposit_l((Word16)ptr[(*bp)++])); move32();
        assert(st_fx->ac_low_fx < (1U << 24));
    }

    st_fx->BER_detect = 0; move16();

    Dyn_Mem_Deluxe_Out();
}

/* o  : Decoded cumulative frequency    */
static __forceinline Word16 ac_decode_fx(Decoder_State_fx *st_fx, /* i/o: Decoder State                   */
                                         Word16            pki)
{
    Dyn_Mem_Deluxe_In(
        UWord16 sgn;
        Word16  val, r;
    );

    st_fx->ac_help_fx = UL_lshr_pos(st_fx->ac_range_fx, 10); move32();
    val               = 0;                                   move16();

    r = add(val, 8);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ari_spec_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r; move16();
    }

    r = add(val, 4);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ari_spec_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r; move16();
    }

    r = add(val, 2);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ari_spec_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r; move16();
    }

    r = add(val, 1);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ari_spec_cumfreq[pki][r]), &sgn);
    IF (sgn == 0)
    {
        val = r; move16();
        IF (sub(val, 15) == 0)
        {
            UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ari_spec_cumfreq[pki][16]), &sgn);
            if (sgn == 0)
            {
                val = 16; move16();
            }
            UL_subNs(st_fx->ac_low_fx, UL_lshl(st_fx->ac_help_fx, 10), &sgn);
            if (sgn == 0)
            {
                st_fx->BER_detect = 1; move16();
            }
        }
    }

    Dyn_Mem_Deluxe_Out();
    return val;
}

/* o  : Decoded cumulative frequency    */
static __forceinline Word16 ac_decode_tns_order(Decoder_State_fx *st_fx, /* i/o: Decoder State                   */
                                                Word16            enable_lpc_weighting)
{
    Dyn_Mem_Deluxe_In(
        UWord16 sgn;
        Word16  val, r;
    );

    st_fx->ac_help_fx = UL_lshr_pos(st_fx->ac_range_fx, 10); move32();
    val               = 0;                                   move16();

    r = add(val, 4);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_order_cumfreq[enable_lpc_weighting][r]), &sgn);
    if (sgn == 0)
    {
        val = r; move16();
    }

    r = add(val, 2);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_order_cumfreq[enable_lpc_weighting][r]), &sgn);
    if (sgn == 0)
    {
        val = r; move16();
    }

    r = add(val, 1);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_order_cumfreq[enable_lpc_weighting][r]), &sgn);
    if (sgn == 0)
    {
        val = r; move16();
    }

    UL_subNs(st_fx->ac_low_fx, UL_lshl(st_fx->ac_help_fx, 10), &sgn);
    if (sgn == 0)
    {
        st_fx->BER_detect = 1; move16();
    }

    Dyn_Mem_Deluxe_Out();
    return val;
}

 /* o  : Decoded cumulative frequency    */
static __forceinline Word16 ac_decode_tns_coef(Decoder_State_fx *st_fx, /* i/o: Decoder State                   */
                                               Word16            pki)
{
    Dyn_Mem_Deluxe_In(
        UWord16 sgn;
        Word16  val, r;
    );

    st_fx->ac_help_fx = UL_lshr_pos(st_fx->ac_range_fx, 10); move32();
    val               = 0;                                   move16();

    r = add(val, 8);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_coef_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r; move16();
    }

    r = add(val, 4);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_coef_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r; move16();
    }

    r = add(val, 2);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_coef_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r; move16();
    }

    r = add(val, 1);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_coef_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r; move16();
        IF (sub(val, 15) == 0)
        {
            UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_coef_cumfreq[pki][16]), &sgn);
            if (sgn == 0)
            {
                val = 16; move16();
            }
            UL_subNs(st_fx->ac_low_fx, UL_lshl(st_fx->ac_help_fx, 10), &sgn);
            if (sgn == 0)
            {
                st_fx->BER_detect = 1; move16();
            }
        }
    }

    Dyn_Mem_Deluxe_Out();
    return val;
}

static __forceinline Word16 ac_dec_update_fx(UWord8 *ptr, Word16 *bp, Word16 *bp_side, Word16 *mask_side, Word16 cur_bin,
                                             Decoder_State_fx *st_fx, /* i/o: Decoder State */
                                             UWord32 cum_freq, /* i  : Cumulative frequency    */
                                             UWord32 sym_freq  /* i  : Symbol frequency        */
)
{
    UWord32 UL_tmp;


    assert(st_fx->ac_help_fx < (1U << 24));
    assert(cum_freq < (1U << 24));

    UL_tmp = UL_Mpy_32_32(cum_freq, st_fx->ac_help_fx);
    assert(UL_tmp < (1U << 24));

    st_fx->ac_low_fx = UL_subNsD(st_fx->ac_low_fx, UL_tmp); move32(); /*0+0*/
    assert(st_fx->ac_low_fx < (1U << 24));

    st_fx->ac_range_fx = UL_Mpy_32_32(st_fx->ac_help_fx, sym_freq); move32();

    assert(st_fx->ac_range_fx < (1U << 24));
    /* updated to 16 from 24 */
    WHILE (st_fx->ac_range_fx < (1U << 16))
    {
        L_sub(0, 0); /* For comparision in while*/

        st_fx->ac_low_fx =
            UL_and(st_fx->ac_low_fx, 0x0000ffFF); /*  make sure upshift doe not lead to more than 24 bits */
        assert(st_fx->ac_low_fx < 1U << 16);

        if (check_pc_bytes(bp, bp_side, mask_side, cur_bin, 1, &st_fx->pc) != 0) return 1;

        /*shift in 8 bits */
        st_fx->ac_low_fx = UL_addNsD(UL_lshl_pos(st_fx->ac_low_fx, 8), UL_deposit_l((Word16)ptr[(*bp)++])); move32();

        assert(st_fx->ac_low_fx < (1U << 24));
        st_fx->ac_range_fx = UL_lshl_pos(st_fx->ac_range_fx, 8); move32();
        assert(st_fx->ac_range_fx < (1U << 24));
    }
    return 0;
}

static __forceinline void pc_init_fx(Word16 n_pc, Word16 numbytes, Word16 be_bp_left, Word16 be_bp_right, Word16 L_spec,
                                     Word16 enc, Word16 bfi, Pc_State_fx *pc /* i/o: Pc State */
)
{
    pc->inv_bin     = add(L_spec, 1);       move16();
    pc->numbytes    = numbytes;             move16();
    pc->c_bp        = 0;                    move16();
    pc->c_bp_side   = 0;                    move16();
    pc->bytes       = shr(add(n_pc, 1),1);  move16();
    pc->b_left      = add(numbytes,1);      move16();
    pc->b_right     = -1;                   move16();
    pc->enc         = enc;                  move16();
    pc->bfi         = bfi;                  move16();
    pc->be_bp_left  = shr(be_bp_left, 3);   move16();
    pc->be_bp_right = shr(be_bp_right, 3);  move16();
    assert(pc->be_bp_right < pc->bytes || pc->bytes == 0);
}

static __forceinline Word16 check_pc_bytes(Word16 *bp, Word16 *bp_side, Word16 *mask_side, Word16 cur_bin, Word16 from_left,
                                           Pc_State_fx *pc /* i/o: Pc State */)
{
    Dyn_Mem_Deluxe_In(
        Word16 bp_local, bp_side_local, offset;
    );

    IF (pc->bytes > 0)
    {
        test();
        IF (from_left == 0 && sub(*mask_side, 1) != 0)
        {
            Dyn_Mem_Deluxe_Out();
            return 0;
        }
        bp_local = *bp;
        bp_side_local = *bp_side;

        IF (from_left != 0)
        {
            if (sub(*mask_side, 1) == 0)
            {
                bp_side_local = add(bp_side_local, 1);
            }
        }
        ELSE
        {
            bp_local = sub(bp_local, 1);
        }

        IF (pc->b_right < 0)
        {
            offset = -1;  move16();
            if (pc->enc == 0)
            {
                offset = add(offset, pc->bytes);
            }

            IF (add(bp_side_local, sub(offset, bp_local)) == pc->bytes)
            {
                pc->b_left = add(bp_local, 1);
                pc->b_right = sub(bp_side_local, 1);
                IF (pc->enc != 0)
                {
                    assert(pc->b_right-pc->b_left+1 == pc->bytes);
                    Dyn_Mem_Deluxe_Out();
                    return 1;
                }
            }
        }

        test();
        IF (pc->enc == 0 && pc->b_right >= 0)
        {
            test();
            IF (from_left != 0 && sub(*bp, pc->b_left) == 0)
            {
                *bp = 0;  move16();
                pc->c_bp = 1;  move16();
            }
            test();
            IF (from_left == 0 && sub(bp_side_local, pc->b_right) == 0)
            {
                *bp_side = sub(pc->bytes, 1);  move16();
                pc->c_bp_side = 1;  move16();
            }
            IF (sub(pc->bfi, 2) == 0)
            {
                test();test();
                IF ((pc->c_bp != 0 && sub(*bp, pc->be_bp_left) >= 0) || (pc->c_bp_side != 0 && sub(*bp_side, pc->be_bp_right) <= 0))
                {
                    pc->inv_bin = cur_bin; move16();
                    Dyn_Mem_Deluxe_Out();
                    return 1;
                }
                ELSE IF((pc->c_bp != 0 && *bp >= 0) || (pc->c_bp_side != 0 && sub(*bp_side, sub(pc->bytes, 1)) <= 0))
                {
                    pc->inv_bin = s_min(pc->inv_bin, cur_bin);
                    Dyn_Mem_Deluxe_Out();
                    return 0;
                }
            }
        }
    }

    Dyn_Mem_Deluxe_Out();
    return 0;
}

