/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static void ac_shift_fl(Encoder_State_fl* st);
static void ac_encode_fl(Encoder_State_fl* st, LC3_INT sym_freq, LC3_INT cum_freq);
static void tns_order_freq_enc(LC3_INT enable_lpc_weighting, LC3_INT order, LC3_INT* symfreq, LC3_INT* cumfreq);
static void tns_coef_freq_enc(LC3_INT k, LC3_INT idx, LC3_INT* symfreq, LC3_INT* cumfreq);
static void ac_freq_fl(LC3_INT pki, LC3_INT s, LC3_INT* symfreq, LC3_INT* cumfreq);
static void ac_finalize_fl(Encoder_State_fl* st);
static void write_uint_forward_fl(Encoder_State_fl* st, LC3_INT val, LC3_INT numbits);
static void ari_enc_init(Encoder_State_fl* st, LC3_UINT8* bytes, LC3_INT* bp_side, LC3_INT* mask_side);
static LC3_INT  sign(LC3_INT x);
static void read_bit_fl(LC3_UINT8* ptr, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT* bit);
static void ac_dec_init_fl(LC3_UINT8* ptr, LC3_INT* bp, Decoder_State_fl* st_fl);
static void tns_order_freq(LC3_INT enable_lpc_weighting, LC3_INT* symfreq, LC3_INT* cumfreq, LC3_INT* numsym);
static void tns_coef_freq(LC3_INT k, LC3_INT* symfreq, LC3_INT* cumfreq, LC3_INT* numsym);
static LC3_INT  ac_decode_fl(Decoder_State_fl* st, LC3_INT* sym_freq, LC3_INT* cum_freq, LC3_INT num_sym, LC3_UINT8* ptr, LC3_INT* bp);
static void ac_freq(LC3_INT pki, LC3_INT* symfreq, LC3_INT* cumfreq, LC3_INT* numsym);
static void findNonZero(LC3_INT* in, LC3_INT* out, LC3_INT len, LC3_INT* outLen);

void ac_dec_init_fl(LC3_UINT8* ptr, LC3_INT* bp, Decoder_State_fl* st_fl)
{
    LC3_INT i = 0;

    st_fl->ac_low_fl = 0;

    st_fl->ac_range_fl = (LC3_UINT32)pow(2, 24) - (LC3_UINT32)1;
    for (i = 0; i < 3; i++) {
        st_fl->ac_low_fl = (st_fl->ac_low_fl << 8) + (LC3_UINT32)ptr[*bp];
        *bp              = *bp + 1;
    }

    st_fl->BER_detect = 0;
}

void tns_order_freq(LC3_INT enable_lpc_weighting, LC3_INT* symfreq, LC3_INT* cumfreq, LC3_INT* numsym)
{
    LC3_INT i = 0, j = 0;

    *numsym = 8;

    j = 0;
    for (i = 1; i < 9; i++) {
        symfreq[j] = ari_tns_order_cf[enable_lpc_weighting][i];
        j++;
    }

    for (i = 0; i < *numsym; i++) {
        symfreq[i] -= ari_tns_order_cf[enable_lpc_weighting][i];
    }

    for (i = 0; i < *numsym; i++) {
        cumfreq[i] = ari_tns_order_cf[enable_lpc_weighting][i];
    }
}

/* Returns val */
LC3_INT ac_decode_fl(Decoder_State_fl* st, LC3_INT* sym_freq, LC3_INT* cum_freq, LC3_INT num_sym, LC3_UINT8* ptr, LC3_INT* bp)
{
    LC3_INT val = 0, tmp = 0;

    tmp = st->ac_range_fl >> 10;

    if (st->ac_low_fl >= (LC3_UINT32)(tmp << 10)) {
        st->BER_detect = 1;
    }

    val = num_sym - 1;

    while (st->ac_low_fl < (LC3_UINT32)(tmp * cum_freq[val])) {
        val--;
    }

    st->ac_low_fl   = st->ac_low_fl - tmp * cum_freq[val];
    st->ac_range_fl = tmp * sym_freq[val];

    while (st->ac_range_fl < pow(2, 16)) {
        st->ac_low_fl   = st->ac_low_fl << 8;
        st->ac_low_fl   = ((LC3_INT)st->ac_low_fl) & ((LC3_INT)(pow(2, 24) - 1));
        st->ac_low_fl   = st->ac_low_fl + ptr[*bp];
        *bp             = *bp + 1;
        st->ac_range_fl = st->ac_range_fl << 8;
    }

    return val;
}

void tns_coef_freq(LC3_INT k, LC3_INT* symfreq, LC3_INT* cumfreq, LC3_INT* numsym)
{
    LC3_INT i = 0, j = 0;

    *numsym = 18 - 1;

    j = 0;
    for (i = 1; i <= *numsym; i++) {
        symfreq[j] = ari_tns_freq_cf[k][i];
        j++;
    }

    for (i = 0; i < *numsym; i++) {
        symfreq[i] -= ari_tns_freq_cf[k][i];
    }

    for (i = 0; i < *numsym; i++) {
        cumfreq[i] = ari_tns_freq_cf[k][i];
    }
}

void ac_freq(LC3_INT pki, LC3_INT* symfreq, LC3_INT* cumfreq, LC3_INT* numsym)
{
    LC3_INT i = 0, j = 0;

    *numsym = 18 - 1;

    j = 0;
    for (i = 1; i <= *numsym; i++) {
        symfreq[j] = ari_spec_cumfreq_fl[pki][i];
        j++;
    }

    for (i = 0; i < *numsym; i++) {
        symfreq[i] -= ari_spec_cumfreq_fl[pki][i];
    }

    for (i = 0; i < *numsym; i++) {
        cumfreq[i] = ari_spec_cumfreq_fl[pki][i];
    }
}

void read_bit_fl(LC3_UINT8* ptr, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT* bit)
{
    *bit = 0;

    if (ptr[*bp_side] & *mask_side) {
        *bit = 1;
    } else {
        *bit = 0;
    }

    if (*mask_side == 128) {
        *mask_side = 1;
        *bp_side   = *bp_side - 1;
    } else {
        *mask_side = *mask_side * 2;
    }
}

void findNonZero(LC3_INT* in, LC3_INT* out, LC3_INT len, LC3_INT* outLen)
{
    LC3_INT i = 0, j = 0;

    for (i = 0; i < len; i++) {
        if (in[i] != 0) {
            out[j] = i;
            j++;
        }
    }

    *outLen = j;
}

void processAriDecoder_fl(LC3_UINT8* bytes, LC3_INT bp_side, LC3_INT mask_side, LC3_INT L_spec, LC3_INT fs_idx, LC3_INT enable_lpc_weighting,
                          LC3_INT tns_numfilters, LC3_INT lsbMode, LC3_INT lastnz, LC3_INT* bfi, LC3_INT* tns_order, LC3_INT fac_ns_idx,
                          LC3_INT gg_idx, uint8_t * resBits, LC3_INT* x, LC3_INT* nf_seed, LC3_INT* tns_idx, LC3_INT* zero_frame, LC3_INT numbytes,
                          LC3_INT* nbits_residual, LC3_INT* residualPresent
#ifdef ENABLE_HR_MODE
                          ,
                          LC3_INT hrmode
#endif
)
{
    Decoder_State_fl st;
    LC3_INT              a = 0, b = 0, t = 0, bp = 0;
    LC3_INT              c = 0;
    LC3_INT              nbits_side = 0;
    LC3_UINT8*           ptr = NULL;
    LC3_INT              n = 0, k = 0, lev = 0;
    LC3_INT              max_lev = 0;
    LC3_INT              sym_freq[MAX_LEN] = {0}, cum_freq[MAX_LEN] = {0}, numsym = 0, bit = 0, lev1 = 0, pki = 0, sym = 0,
        save_lev[MAX_LEN] = {0}, idx_len = 0, total_bits = 0, nbits_ari = 0, idx[MAX_LEN] = {0}, rateFlag = 0;

    total_bits = 8 * numbytes;

    /* Rate flag */
    if (fs_idx != 5)
    {
        if (total_bits > (160 + fs_idx * 160)) {
            rateFlag = 512;
        }
    }

    /* Init */
    c  = 0;
    t  = 0;
    bp = 0;

    ptr = bytes;

    /* Start Decoding */
    ac_dec_init_fl(ptr, &bp, &st);

    /* Decode TNS data */
    for (n = 0; n < tns_numfilters; n++) {
        if (tns_order[n] > 0) {
            tns_order_freq(enable_lpc_weighting, sym_freq, cum_freq, &numsym);
            tns_order[n] = ac_decode_fl(&st, sym_freq, cum_freq, numsym, ptr, &bp);
            tns_order[n] = tns_order[n] + 1;

            for (k = 0; k < tns_order[n]; k++) {
                tns_coef_freq(k, sym_freq, cum_freq, &numsym);
                tns_idx[n * 8 + k] = ac_decode_fl(&st, sym_freq, cum_freq, numsym, ptr, &bp);
            }
        }
    }

    if (st.BER_detect > 0) {
        *bfi = 1;
        
        return;
    }

    /* Spectral data */
    for (k = 0; k < lastnz; k = k + 2) {
        /* Context */
        t = c + rateFlag;

        if (k > L_spec / 2) {
            t = t + 256;
        }

        /* Decode amplitude */
        x[k]     = 0;
        x[k + 1] = 0;

#ifdef ENABLE_HR_MODE
        if (hrmode == 1) {
            max_lev = 13 + 8;
        } else {
            max_lev = 13;
        }
#else
        max_lev = 13;
#endif

        for (lev = 0; lev <= max_lev; lev++) {
            lev1 = MIN(lev, 3);
            pki  = ari_spec_lookup_fl[t + lev1 * 1024];

            ac_freq(pki, sym_freq, cum_freq, &numsym);
            sym = ac_decode_fl(&st, sym_freq, cum_freq, numsym, ptr, &bp);

            if (sym < 16) {
                break;
            }

            if (lsbMode == 0 || lev > 0) {
                read_bit_fl(ptr, &mask_side, &bp_side, &bit);

                x[k] = x[k] + (bit << lev);

                read_bit_fl(ptr, &mask_side, &bp_side, &bit);

                x[k + 1] = x[k + 1] + (bit << lev);
            }
        }

        if (lsbMode == 1) {
            save_lev[k] = lev;
        }

        a = sym & 3;
        b = sym >> 2;

        x[k]     = x[k] + (a << lev);
        x[k + 1] = x[k + 1] + (b << lev);

        /* Decode signs */
        if (x[k] > 0) {
            read_bit_fl(ptr, &mask_side, &bp_side, &bit);

            if (bit == 1) {
                x[k] = -x[k];
            }
        }

        if (x[k + 1] > 0) {
            read_bit_fl(ptr, &mask_side, &bp_side, &bit);
            if (bit == 1) {
                x[k + 1] = -x[k + 1];
            }
        }

        /* Context */
        lev1 = MIN(lev, 3);
        if (lev1 <= 1) {
            t = 1 + (a + b) * (lev1 + 1);
        } else {
            t = 12 + lev1;
        }

        c = (c & 15) * 16 + t;

        if ((bp - bp_side) > 3 || st.BER_detect > 0) {
            *bfi = 1;
            return;
        }
    }

    /* Residual bits */
    nbits_side      = total_bits - (8 * (bp_side + 1) + 8 - LC3_LOG2(mask_side));
    nbits_ari       = (bp + 1 - 3) * 8 + 25 - floor(LC3_LOG2(st.ac_range_fl));
    *nbits_residual = total_bits - (nbits_side + nbits_ari);

    if (*nbits_residual < 0) {
        *bfi = 1;
        return;
    }

    if (lsbMode == 0) {
        findNonZero(x, idx, L_spec, &idx_len);
#ifdef ENABLE_HR_MODE
        if (hrmode)
        {
        	idx_len *= EXT_RES_ITER_MAX;
        }
#endif
        *nbits_residual  = MIN(*nbits_residual, idx_len);
        *residualPresent = 1;

        memset(resBits, 0, MAX_RESBITS_LEN);

        for (k = 0; k < *nbits_residual; k++) {
        	LC3_INT tmp;
        	read_bit_fl(ptr, &mask_side, &bp_side, &tmp);
        	resBits[k >> 3] |= tmp << (k & 7);
        }
    } else {
        for (k = 0; k < lastnz; k = k + 2) {
            if (save_lev[k] > 0) {
                if (*nbits_residual == 0) {
                    break;
                }

                read_bit_fl(ptr, &mask_side, &bp_side, &bit);
                *nbits_residual = *nbits_residual - 1;

                if (bit == 1) {
                    if (x[k] > 0) {
                        x[k] = x[k] + 1;
                    } else if (x[k] < 0) {
                        x[k] = x[k] - 1;
                    } else {
                        if (*nbits_residual == 0) {
                            break;
                        }

                        read_bit_fl(ptr, &mask_side, &bp_side, &bit);
                        *nbits_residual = *nbits_residual - 1;

                        if (bit == 0) {
                            x[k] = 1;
                        } else {
                            x[k] = -1;
                        }
                    }
                }

                if (*nbits_residual == 0) {
                    break;
                }

                read_bit_fl(ptr, &mask_side, &bp_side, &bit);
                *nbits_residual = *nbits_residual - 1;

                if (bit == 1) {
                    if (x[k + 1] > 0) {
                        x[k + 1] = x[k + 1] + 1;
                    } else if (x[k + 1] < 0) {
                        x[k + 1] = x[k + 1] - 1;
                    } else {
                        if (*nbits_residual == 0) {
                            break;
                        }

                        read_bit_fl(ptr, &mask_side, &bp_side, &bit);
                        *nbits_residual = *nbits_residual - 1;

                        if (bit == 0) {
                            x[k + 1] = 1;
                        } else {
                            x[k + 1] = -1;
                        }
                    }
                }
            }
        }
    }

    /* Noise-filling seed */
    *nf_seed = 0;

    for (k = 0; k < L_spec; k++) {
        *nf_seed = *nf_seed + abs(x[k]) * k;
    }

    *nf_seed = *nf_seed & 65535;

    if (*nf_seed >= 32768) {
        *nf_seed = *nf_seed - 65536;
    }

    /* Zero frame flag */
    if (lastnz == 2 && x[0] == 0 && x[1] == 0 && gg_idx == 0 && fac_ns_idx == 7) {
        *zero_frame = 1;
    } else {
        *zero_frame = 0;
    }
}

void ac_encode_fl(Encoder_State_fl* st, LC3_INT sym_freq, LC3_INT cum_freq)
{
    LC3_INT r = 0;

    r       = st->range >> 10;
    st->low = st->low + r * cum_freq;

    if ((st->low >> 24) == 1) {
        st->carry = 1;
    }

    st->low   = (st->low) & ((LC3_INT)pow(2, 24) - 1);
    st->range = r * sym_freq;

    while (st->range < (LC3_INT)pow(2, 16)) {
        st->range = st->range << 8;
        ac_shift_fl(st);
    }
}

void ac_shift_fl(Encoder_State_fl* st)
{
    if (st->low < 16711680 || st->carry == 1) {
        if (st->cache >= 0) {
            st->ptr[st->bp] = st->cache + st->carry;
            st->bp          = st->bp + 1;
        }

        while (st->carry_count > 0) {
            st->ptr[st->bp] = (st->carry + 255) & 255;
            st->bp          = st->bp + 1;
            st->carry_count = st->carry_count - 1;
        }

        st->cache = st->low >> 16;
        st->carry = 0;
    } else {
        st->carry_count = st->carry_count + 1;
    }

    st->low = st->low << 8;
    st->low = (st->low) & ((LC3_INT)pow(2, 24) - 1);
}

void tns_order_freq_enc(LC3_INT enable_lpc_weighting, LC3_INT order, LC3_INT* symfreq, LC3_INT* cumfreq)
{
    *symfreq = tns_freq_cf[enable_lpc_weighting][order] - tns_freq_cf[enable_lpc_weighting][order - 1];
    *cumfreq = tns_freq_cf[enable_lpc_weighting][order - 1];
}

void tns_coef_freq_enc(LC3_INT k, LC3_INT idx, LC3_INT* symfreq, LC3_INT* cumfreq)
{
    *symfreq = tns_cf[k][idx + 1] - tns_cf[k][idx];
    *cumfreq = tns_cf[k][idx];
}

void ac_freq_fl(LC3_INT pki, LC3_INT s, LC3_INT* symfreq, LC3_INT* cumfreq)
{
    *symfreq = ari_spec_cumfreq_fl[pki][s + 1] - ari_spec_cumfreq_fl[pki][s];
    *cumfreq = ari_spec_cumfreq_fl[pki][s];
}

void ac_finalize_fl(Encoder_State_fl* st)
{
    LC3_INT bits = 0, mask = 0, val = 0, over1 = 0, high = 0, over2 = 0, c = 0, b = 0;

    bits  = 24 - floor(LC3_LOG2(st->range));
    mask  = ((LC3_INT)pow(2, 24) - 1) >> bits;
    val   = st->low + mask;
    over1 = val >> 24;

    val   = (val) & ((LC3_INT)pow(2, 24) - 1);
    high  = st->low + st->range;
    over2 = high >> 24;
    high  = high & ((LC3_INT)pow(2, 24) - 1);
    val   = val & (((LC3_INT)pow(2, 24) - 1) - mask);

    if (over1 == over2) {
        if (val + mask >= high) {
            bits = bits + 1;
            mask = mask >> 1;
            val  = ((st->low + mask) & ((LC3_INT)pow(2, 24) - 1)) & (((LC3_INT)pow(2, 24) - 1) - mask);
        }

        if (val < st->low) {
            st->carry = 1;
        }
    }

    st->low = val;

    b = bits;

    if (bits > 8) {
        for (; b >= 1; b = b - 8) {
            ac_shift_fl(st);
        }
    } else {
        ac_shift_fl(st);
    }

    bits = b;
    if (bits < 0) {
        bits += 8;
    }

    if (st->carry_count > 0) {
        st->ptr[st->bp] = st->cache;
        st->bp          = st->bp + 1;

        for (c = st->carry_count; c >= 2; c--) {
            st->ptr[st->bp] = 255;
            st->bp          = st->bp + 1;
        }

        write_uint_forward_fl(st, 255 << (bits - 8), bits);
    } else {
        write_uint_forward_fl(st, st->cache, bits);
    }
}

void write_uint_forward_fl(Encoder_State_fl* st, LC3_INT val, LC3_INT numbits)
{
    LC3_INT k = 0, bit = 0, mask = 128;

    for (k = 0; k < numbits; k++) {
        bit = val & mask;

        if (bit == 0) {
            st->ptr[st->bp] = st->ptr[st->bp] & (255 - mask);
        } else {
            st->ptr[st->bp] = st->ptr[st->bp] | mask;
        }

        mask = mask >> 1;
    }
}

void ari_enc_init(Encoder_State_fl* st, LC3_UINT8* bytes, LC3_INT* bp_side, LC3_INT* mask_side)
{
    st->ptr         = bytes;
    st->bp_side     = bp_side;
    st->mask_side   = mask_side;
    st->bp          = 0;
    st->low         = 0;
    st->range       = (LC3_INT)pow(2, 24) - 1;
    st->cache       = -1;
    st->carry       = 0;
    st->carry_count = 0;
}

LC3_INT sign(LC3_INT x)
{
    if (x > 0)
        return 1;

    if (x < 0)
        return -1;

    return 0;
}

void processAriEncoder_fl(LC3_UINT8* bytes, LC3_INT bp_side, LC3_INT mask_side, LC3_INT* x, LC3_INT* tns_order, LC3_INT tns_numfilters,
                          LC3_INT* tns_idx, LC3_INT lastnz, LC3_INT* codingdata, uint8_t* res_bits, LC3_INT resBitsLen, LC3_INT lsbMode,
                          LC3_INT nbbits, LC3_INT enable_lpc_weighting)
{
    LC3_INT              total_bits = 0, cumfreq = 0, symfreq = 0, k = 0, i = 0, j = 0, lev = 0, lev1 = 0;
    LC3_INT              bit1 = 0, bit2 = 0, lsb1 = 0, lsb2 = 0, a = 0, b = 0, bit = 0, pki = 0, nbits_side = 0;
    LC3_INT              nbits_residual_enc = 0, nbits_ari = 0, lsbs[MAX_LEN] = {0}, lsbsLen = 0;
    Encoder_State_fl st;

    ari_enc_init(&st, bytes, &bp_side, &mask_side);

    total_bits = nbbits;

    /* TNS data  */
    for (i = 0; i < tns_numfilters; i++) {
        if (tns_order[i] > 0) {
            tns_order_freq_enc(enable_lpc_weighting, tns_order[i], &symfreq, &cumfreq);
            ac_encode_fl(&st, symfreq, cumfreq);

            for (j = 0; j < tns_order[i]; j++) {
                tns_coef_freq_enc(j, tns_idx[i * 8 + j], &symfreq, &cumfreq);
                ac_encode_fl(&st, symfreq, cumfreq);
            }
        }
    }

    /* Spectral data */
    for (k = 0; k < lastnz; k = k + 2) {
        for (lev = 0; lev < codingdata[1]; lev++) {
            lev1 = MIN(lev, 3);
            pki  = ari_spec_lookup_fl[codingdata[0] + lev1 * 1024];
            ac_freq_fl(pki, 16, &symfreq, &cumfreq);

            ac_encode_fl(&st, symfreq, cumfreq);
            bit1 = (abs(x[k]) >> lev) & 1;
            bit2 = (abs(x[k + 1]) >> lev) & 1;

            if (lsbMode == 1 && lev == 0) {
                lsb1 = bit1;
                lsb2 = bit2;
            } else {
                write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, bit1);
                write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, bit2);
            }
        }

        lev1 = MIN(MAX(codingdata[1], 0), 3);
        pki  = ari_spec_lookup_fl[codingdata[0] + lev1 * 1024];

        ac_freq_fl(pki, codingdata[2], &symfreq, &cumfreq);
        ac_encode_fl(&st, symfreq, cumfreq);

        a = abs(x[k]);
        b = abs(x[k + 1]);

        if (lsbMode == 1 && codingdata[1] > 0) {
            a             = a >> 1;
            lsbs[lsbsLen] = lsb1;
            lsbsLen++;

            if (a == 0 && x[k] != 0) {
                bit           = MAX(0, -sign(x[k]));
                lsbs[lsbsLen] = bit;
                lsbsLen++;
            }

            b             = b >> 1;
            lsbs[lsbsLen] = lsb2;
            lsbsLen++;

            if (b == 0 && x[k + 1] != 0) {
                bit           = MAX(0, -sign(x[k + 1]));
                lsbs[lsbsLen] = bit;
                lsbsLen++;
            }
        }

        if (a != 0) {
            bit = MAX(0, -sign(x[k]));
            write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, bit);
        }

        if (b != 0) {
            bit = MAX(0, -sign(x[k + 1]));
            write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, bit);
        }

        codingdata += 3;
    }

    /* Residual bits */
    nbits_side = total_bits - (8 * (*(st.bp_side) + 1) + 8 - LC3_LOG2(*(st.mask_side)));
    nbits_ari  = (st.bp + 1) * 8 + 25 - floor(LC3_LOG2(st.range));

    if (st.cache >= 0) {
        nbits_ari = nbits_ari + 8;
    }

    if (st.carry_count > 0) {
        nbits_ari = nbits_ari + st.carry_count * 8;
    }

    nbits_residual_enc = total_bits - (nbits_side + nbits_ari);
    
    assert(nbits_residual_enc >= 0);

    if (lsbMode == 0) {
        nbits_residual_enc = MIN(nbits_residual_enc, resBitsLen);
        for (k = 0; k < nbits_residual_enc; k++) {
        	if (res_bits[k >> 3] & (1 << (k & 7)))
        	{
        		write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, 1);
        	}
        	else
        	{
        		write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, 0);
        	}
        }
    } else {
        nbits_residual_enc = MIN(nbits_residual_enc, lsbsLen);

        for (k = 0; k < nbits_residual_enc; k++) {
            write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, lsbs[k]);
        }
    }

    ac_finalize_fl(&st);
}
