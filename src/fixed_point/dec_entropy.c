/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "functions.h"


static Word16 read_indice(UWord8 *ptr, Word16 *bp, Word16 *mask, Word16 numbits);

static Word16 ac_dec_split_st2VQ_CW(                     /* local BER flag */
                                    const Word32 L_cwRx, /* max 25 bits */
                                    const Word32 L_szA, const Word32 L_szB, Word32 *L_cwA, Word32 *L_cwB,
                                    Word16 *submodeLSB);

void processDecoderEntropy_fx(UWord8 *bytes, Word16 *bp_side, Word16 *mask_side, Word16 nbbits,
                              Word16 L_spec, Word16 fs_idx, Word16 BW_cutoff_bits, Word16 *tns_numfilters,
                              Word16 *lsbMode, Word16 *lastnz, Word16 *bfi, Word16 *tns_order, Word16 *fac_ns_idx,
                              Word16 *gg_idx, Word16 *BW_cutoff_idx, Word16 *ltpf_idx, Word32 *L_scf_idx,
                              Word16 frame_dms)
{
    Dyn_Mem_Deluxe_In(
        Word16  L, submodeLSB;
        Word32  tmp32, tmp32lim;
        Word16  gain_e, gain, submodeMSB, BER_detect;
        Counter n;
        UWord8 *ptr;
    );

    ptr        = bytes;
    *bp_side   = shr_pos(sub(nbbits, 1), 3);
    *mask_side = shl(1, sub(8, sub(nbbits, shl_pos(*bp_side, 3))));

    /* Cutoff-detection */
    IF (BW_cutoff_bits > 0)
    {
        *BW_cutoff_idx = read_indice(ptr, bp_side, mask_side, BW_cutoff_bits);
        /* check for bitflips */
        IF (sub(fs_idx, *BW_cutoff_idx) < 0)
        {
            *BW_cutoff_idx = fs_idx;
            *bfi           = 1;  move16();
            Dyn_Mem_Deluxe_Out();
            return;
        }
    }
    ELSE
    {
        *BW_cutoff_idx = 0;
    }

    /* Number of TNS filters */
    IF (sub(*BW_cutoff_idx, 3) >= 0 && frame_dms >= 50)
    {
        *tns_numfilters = 2;  move16();
    }
    ELSE
    {
        *tns_numfilters = 1;  move16();
    }

    /* Decode number of ntuples */
    L       = sub(14, norm_s(negate(L_spec)));
    n       = read_indice(ptr, bp_side, mask_side, L);
    n       = add(n, 1);
    *lastnz = shl_pos(n, 1);
    IF (sub(*lastnz, L_spec) > 0)
    {
        *bfi = 1;  move16();
        Dyn_Mem_Deluxe_Out();
        return;
    }

    /* Mode bit */
    *lsbMode = read_bit(ptr, bp_side, mask_side);

    /* Decode global-gain */
    *gg_idx = read_indice(ptr, bp_side, mask_side, 8);  move16();
    tmp32  = L_shl_pos(L_mult0(*gg_idx, 0x797D), 7);  /* 6Q25; 0x797D -> log2(10)/28 (Q18) */
    gain_e = add(extract_l(L_shr_pos(tmp32, 25)), 1); /* get exponent */
    gain   = round_fx(BASOP_Util_InvLog2(L_or(tmp32, 0xFE000000))); 
    assert(gain >= 0); /* JSv, check if shr_pos(gain,1)  is more appropriate) */
    gain   = shr_r(gain, 1);
    gain_e = add(gain_e, 1);

    /* Decode TNS on/off flag */
    tns_order[1] = 0; move16(); /* fix problem with uninitialized memory */
    FOR (n = 0; n < *tns_numfilters; n++)
    {
        tns_order[n] = read_bit(ptr, bp_side, mask_side);  move16();
    }

    /* LTPF on/off */
    ltpf_idx[0] = read_indice(ptr, bp_side, mask_side, 1);  move16();

    /* Decode SNS VQ parameters - 1st stage (10 bits) */
    L_scf_idx[0] = L_deposit_l(read_indice(ptr, bp_side, mask_side, 5)); /* stage1 LF  5  bits */
    L_scf_idx[1] = L_deposit_l(read_indice(ptr, bp_side, mask_side, 5)); /* stage1 HF  5 bits  */

    /* Decode SNS VQ parameters - 2nd stage side-info (3-4 bits) */
    submodeMSB   = read_bit(ptr, bp_side, mask_side); /* submodeMSB 1 bit */
    L_scf_idx[2] = L_deposit_l(shl_pos(submodeMSB, 1));
    ASSERT(sns_gainMSBbits[L_scf_idx[2]] > 0);
    L_scf_idx[3] = L_deposit_l(
        read_indice(ptr, bp_side, mask_side, sns_gainMSBbits[L_scf_idx[2]])); /* gains or gain MSBs  1-2 bits  */
    L_scf_idx[4] = read_bit(ptr, bp_side, mask_side);                         /*  shape LS 1 bit */

    /* Decode SNS VQ parameters - 2nd stage data (24-25 bits) */
    IF (submodeMSB == 0)
    { /* shape_j = 0, or 1  */
        /* regular mode A,B indexes integer multiplexed, total 24.x bits  MPVQ codeword section A and  codeword for
         * section B */
        /* regular mode  mode shape  index   total  24.9999 bits    MPVQ codeword  */
        tmp32 = L_deposit_l(read_indice(ptr, bp_side, mask_side, 13));
        tmp32 = L_or(tmp32, L_shl_pos(read_indice(ptr, bp_side, mask_side, 12), 13));  move16(); /*for ber state   */
        BER_detect =
            ac_dec_split_st2VQ_CW(       /* local BER flag */
                                  tmp32, /* L_cwRx  max 25 bits */
                                  sns_MPVQ_Sz[0][0], UL_addNsD(sns_MPVQ_Sz[0][1], sns_MPVQ_Sz[1][1]), /* 12+2 = 14 */
                                  (&L_scf_idx[5]),                                                    /* shape A */
                                  (&L_scf_idx[6]), /* shape B or  gain LSB */
                                  &submodeLSB      /* total submode update below  */
            );
        IF (submodeLSB != 0)
        { /* add gainLSB bit */
            L_scf_idx[3] = L_add(L_shl_pos(L_scf_idx[3], 1), L_scf_idx[6]);
            L_scf_idx[6] = -2L;
        }
    }
    ELSE
    { /* shape_j = 2 or 3  */
        ASSERT(submodeMSB == 1);
        /* outlier mode shape  index   total  23.8536 +  19.5637 (19.5637 < (log2(2.^24 -2.^23.8537))    bits    MPVQ
         * codeword  */
        tmp32        = L_deposit_l(read_indice(ptr, bp_side, mask_side, 12));
        tmp32        = L_or(tmp32, L_shl_pos(read_indice(ptr, bp_side, mask_side, 12), 12));
        L_scf_idx[5] = tmp32;  move32(); /*shape outl_near or outl_far */
        submodeLSB = 0;  move16();
        BER_detect = 0;  move16();
        tmp32lim = L_add(sns_MPVQ_Sz[2][0], L_shl_pos(sns_MPVQ_Sz[3][0], 1));
        IF (L_sub(tmp32, tmp32lim) >= 0)
        {
            BER_detect = 1;  move16();
        }
        ELSE
        {
            tmp32 = L_sub(tmp32, sns_MPVQ_Sz[2][0]); /*  a potential high index is computed */
            IF (tmp32 >= 0)
            {
                submodeLSB = 1;  move16();
                ASSERT(tmp32 >= 0 && tmp32 < (Word32)(2 * sns_MPVQ_Sz[3][0]));
                L_scf_idx[3] = L_add(L_shl_pos(L_scf_idx[3], 1), L_and(tmp32, 0x1)); /* add LSB_gain bit to gain MSBs */
                L_scf_idx[5] = L_shr_pos(tmp32, 1); /* MPVQ index with offset and gainLSB removed */
                L_scf_idx[6] = -2L;  move32();
            }
            ELSE
            {
                L_scf_idx[6] = -1L;  move32();
            }
        }
    }
    L_scf_idx[2] =
        L_add(L_scf_idx[2], L_deposit_l(submodeLSB)); /* decoder internal signal shape_j = submode 0..3 to VQ */

    IF (BER_detect > 0)
    {
        *bfi = 1;  move16();
        Dyn_Mem_Deluxe_Out();
        return;
    }

    /* LTPF data */
    IF (ltpf_idx[0] != 0)
    {
        ltpf_idx[1] = read_indice(ptr, bp_side, mask_side, 1);  move16();
        ltpf_idx[2] = read_indice(ptr, bp_side, mask_side, 9);  move16();
    }
    ELSE
    {
        ltpf_idx[1] = 0;  move16();
        ltpf_idx[2] = 0;  move16();
    }

    /* Decode noise-fac */
    *fac_ns_idx = read_indice(ptr, bp_side, mask_side, 3);  move16();

    Dyn_Mem_Deluxe_Out();
}

#ifdef ENABLE_PADDING
int paddingDec_fx(UWord8 *bytes, Word16 nbbits, Word16 L_spec, Word16 BW_cutoff_bits, Word16 ep_enabled,
                   Word16 *total_padding, Word16 *np_zero)
{
    Word16 lastnz_threshold;
    Word16 padding_len_bits, padding_len;

    Word16 bp_side;
    Word16 nbbytes = shr(nbbits,3);

    Word16  mask_side;
    UWord8 *ptr = bytes;

    Word16 lastnz;
    Word16 nbits = sub(14, norm_s(negate(L_spec)));
    *np_zero     = 0;

    *total_padding = 0;

    bp_side   = shr_pos(sub(nbbits, 1), 3);
    mask_side = shl(1, sub(8, sub(nbbits, shl_pos(bp_side, 3))));

    test();
    IF (sub(bp_side, 19) < 0 || sub(bp_side, LC3_MAX_BYTES ) >= 0) {
        return 1;
    }

    ptr = bytes;

    IF (BW_cutoff_bits > 0)
    {
        read_indice(ptr, &bp_side, &mask_side, BW_cutoff_bits);
        move16();
    }

    lastnz = read_indice(ptr, &bp_side, &mask_side, nbits);
    move16();

    lastnz_threshold = sub(shl(1, nbits), 2);

    WHILE (lastnz == lastnz_threshold)
    {
        padding_len_bits = sub(sub(12, nbits), BW_cutoff_bits);

        /*Read padding length*/
        padding_len = read_indice(ptr, &bp_side, &mask_side, padding_len_bits);
        move16();

        /* Read 4 reserved bits */
        read_indice(ptr, &bp_side, &mask_side, 4);
        move16();

        IF (ep_enabled == 0)
        {
            /* Discard padding length bytes */
            bp_side        = sub(bp_side, padding_len);
            *total_padding = add(add(*total_padding, padding_len), 2); move16();
        }
        ELSE
        {
            *total_padding = add(*total_padding, 2); move16();
            *np_zero       = add(*np_zero, padding_len); move16();
        }
        
        /* test if we have less than 20 bytes left; if so frame is broken */
        IF (sub(sub(nbbytes,add(*total_padding,*np_zero)),20) < 0) {
            return 1;
        }

        /* Read bandwidth bits */
        IF (BW_cutoff_bits > 0)
        {
            read_indice(ptr, &bp_side, &mask_side, BW_cutoff_bits);
            move16();
        }

        lastnz = read_indice(ptr, &bp_side, &mask_side, nbits);
        move16();
    }

    IF (ep_enabled != 0)
    {
        *total_padding = add(*total_padding, *np_zero); move16();
    }
    return 0;
}
#endif

static __forceinline Word16 read_indice(UWord8 *ptr, Word16 *bp, Word16 *mask, Word16 numbits)
{
    Dyn_Mem_Deluxe_In(
        Word16  indice, bit;
        Counter i;
    );

    indice = read_bit(ptr, bp, mask);

    FOR (i = 1; i < numbits; i++)
    {
        bit    = read_bit(ptr, bp, mask);
        indice = add(indice, lshl_pos(bit, i));
    }

    Dyn_Mem_Deluxe_Out();
    return indice;
}

static __forceinline Word16 ac_dec_split_st2VQ_CW(                     /* local BER flag */
                                                  const Word32 L_cwRx, /* max 25 bits */
                                                  const Word32 L_szA, const Word32 L_szB, Word32 *L_cwA, Word32 *L_cwB,
                                                  Word16 *submodeLSB)
{
    /* demultiplex:  L_cwRx =   L_cwB(21.z bits) * L_szA(3.y bits)   + L_cwA(21.x bits)); */
    Word16  start, fin, ind;
    Word32  L_tmp, L_max_size;
    Counter i;

    L_max_size = (Word32)UL_Mpy_32_32((UWord32)L_szB, (UWord32)L_szA); /*  may be tabled  */

    /* section B  ind larger than 13  out of the possible  14 =   0..13  */
    IF (L_sub(L_cwRx, L_max_size) >= 0)
    {
        *L_cwA      = L_deposit_l(0);
        *L_cwB      = L_deposit_l(0);
        *submodeLSB = 0;  move16();
        return (Word16)1; /* set berFlag and exit */
    }

    /*initial binary split of cw,  select top or low half */
    start = 0;  move16();

    ASSERT((L_szB & 0x1L) == 0); /* this middle split only works if  L_szB is even  */
    if (L_sub(L_cwRx, L_shr_pos(L_max_size, 1)) >= 0)
    {
        start = L_shr_pos(L_szB, 1); /* top half start index */
    }

    /*linear loop over a low  or a  high section */
    ind = start;  move16();
    L_tmp = L_negate(L_cwRx); /* search from negative side */

    L_tmp = L_add(L_tmp, (Word32)UL_Mpy_32_32(UL_deposit_l((UWord16)start), (UWord32)L_szA));
    /* start is 0 or 7 */ /*non-fractional mult is   (int)start * L_szA */

    /* a short linear run  over  ceil(szB/2) =  7   values  */

    fin = add(start, shr_pos(L_szB, 1));
    FOR (i = start; i < fin; i++)
    {
        ind   = add(ind, 1);
        L_tmp = L_add(L_tmp, L_szA);
        if (L_tmp > 0)
        {
            ind = sub(ind, 1); /* passed criteria point, keep index    */
        }
    }

    *L_cwB = L_deposit_l(ind);
    *L_cwA = L_sub(L_cwRx, (Word32)UL_Mpy_32_32(UL_deposit_l((UWord16)ind),
                                                (UWord32)L_szA)); /* non-fractional mult;   (int)ind * L_szA */

    ASSERT(*L_cwA >= 0 && *L_cwA < L_szA);
    ASSERT(*L_cwB >= 0 && *L_cwB < L_szB);

    *submodeLSB = 0;
    *L_cwB      = L_sub(*L_cwB, 2);
    if (*L_cwB < 0)
    {
        *submodeLSB = 1;  move16();
    }
    *L_cwB = L_mac0(*L_cwB, 2, *submodeLSB); /* add back gain ind if needed */

    return 0; /* no BER */
}
