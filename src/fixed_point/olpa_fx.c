/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

/*************************************************************************/


void process_olpa_fx(Word16 *mem_s6k4_exp, Word16 mem_s12k8[], Word16 mem_s6k4[], Word16 *pitch, Word16 *s12k8,
                     Word16 len, Word16 *normcorr, Word16 *mem_pitch, Word16 s12k8_exp, Word8 *scratchBuffer)
{
    Word32  sum, sum0, sum1, sum2, prod, inv;
    Word16  shift, s6k4_exp, prod_exp, min_pitch, max_pitch;
    Word16  scale0, scale1, scale2, pitch2, normcorr2, len2;
    Word32  max32;
    Word32 *ac;
    Word16 *s6k4;
    Counter n;

    Counter m;
    Word32  L_tmp, L_tmp2;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("process_olpa_fx", sizeof(struct {
                   Word32  sum, sum0, sum1, sum2, prod, inv;
                   Word16  shift, s6k4_exp, prod_exp, min_pitch, max_pitch;
                   Word16  scale0, scale1, scale2, pitch2, normcorr2, len2;
                   Word32  max32;
                   Word32 *ac;
                   Word16 *s6k4;
                   Counter n;
                   Word32  sums[3];
                   Counter m;
                   Word32  L_tmp, L_tmp2;
               }));
#endif

    /* Buffer alignment */
    ac = (Word32 *)scratchAlign(scratchBuffer, 0); /* Size = 4 * RANGE_PITCH_6K4 = 392 bytes */

    /* Downsample input signal by a factor of 2 (12.8kHz -> 6.4kHz) */
    s6k4    = mem_s6k4 + MAX_PITCH_6K4;
    sum     = L_mac(L_mac(L_mult(mem_s12k8[0], 4053), mem_s12k8[1], 7712), mem_s12k8[2], 9239);
    *s6k4++ = round_fx(L_mac(L_mac(sum, s12k8[0], 7712), s12k8[1], 4053)); move16();
    sum     = L_mac(L_mac(L_mult(mem_s12k8[2], 4053), s12k8[0], 7712), s12k8[1], 9239);
    *s6k4++ = round_fx(L_mac(L_mac(sum, s12k8[2], 7712), s12k8[3], 4053)); move16();

    FOR (n = 5; n < len; n += 2)
    {
        sum     = L_mac(L_mac(L_mult(s12k8[n - 4], 4053), s12k8[n - 3], 7712), s12k8[n - 2], 9239);
        *s6k4++ = round_fx(L_mac(L_mac(sum, s12k8[n - 1], 7712), s12k8[n], 4053)); move16();
    }

    mem_s12k8[0] = s12k8[len - 3]; move16();
    mem_s12k8[1] = s12k8[len - 2]; move16();
    mem_s12k8[2] = s12k8[len - 1]; move16();
    len2         = shr(len, 1);

    /* Scale downsampled signal */
    s6k4          = mem_s6k4 + MAX_PITCH_6K4;
    scale0        = sub(getScaleFactor16_0(mem_s6k4, MAX_PITCH_6K4), 3);
    *mem_s6k4_exp = sub(*mem_s6k4_exp, scale0); move16();
    scale1        = sub(getScaleFactor16_0(s6k4, len2), 3);
    s6k4_exp      = sub(s12k8_exp, scale1);
    scale2        = sub(*mem_s6k4_exp, s6k4_exp);
    IF (scale2 > 0)
    {
        Scale_sig(s6k4, len2, sub(scale1, scale2));
        shift    = scale0;        move16();
        s6k4_exp = *mem_s6k4_exp; move16();
    }
    ELSE
    {
        Scale_sig(s6k4, len2, scale1);
        shift         = add(scale0, scale2);
        *mem_s6k4_exp = s6k4_exp; move16();
    }
    Scale_sig(mem_s6k4, MAX_PITCH_6K4, shift);

    /* Compute autocorrelation */
    FOR (n = MIN_PITCH_6K4; n <= MAX_PITCH_6K4; n++)
    {
        sum = L_mult0(s6k4[0], s6k4[0 - n]);
        FOR (m = 1; m < len2; m++)
        {
            sum = L_mac0(sum, s6k4[m], s6k4[m - n]);
        }
        ac[n - MIN_PITCH_6K4] = sum; move32();
    }

    /* Weight autocorrelation and find maximum */
    max32  = Mpy_32_16(ac[0], olpa_ac_weighting[0]); move32();
    *pitch = MIN_PITCH_6K4;                          move16();
    FOR (n = MIN_PITCH_6K4 + 1; n <= MAX_PITCH_6K4; n++)
    {
        L_tmp  = Mpy_32_16(ac[n - MIN_PITCH_6K4], olpa_ac_weighting[n - MIN_PITCH_6K4]);
        L_tmp2 = L_sub_sat(L_tmp, max32);
        if (L_tmp2 > 0)
        {
            *pitch = n; move16();
        }
        max32 = L_max(L_tmp, max32);
    }

/* Compute normalized correlation */
    sum0 = L_mult0(s6k4[0], s6k4[0 - *pitch]);
    sum1 = L_mac0(1, s6k4[0 - *pitch], s6k4[0 - *pitch]);
    sum2 = L_mac0(1, s6k4[0], s6k4[0]);
    for (m = 1; m < len2; m++)
    {
        sum0 = L_mac0(sum0, s6k4[m], s6k4[m - *pitch]);
        sum1 = L_mac0(sum1, s6k4[m - *pitch], s6k4[m - *pitch]);
        sum2 = L_mac0(sum2, s6k4[m], s6k4[m]);
    }
    scale1   = norm_l(sum1);
    scale2   = norm_l(sum2);
    sum1     = L_shl_pos(sum1, scale1);
    sum2     = L_shl_pos(sum2, scale2);
    prod     = Mpy_32_32(sum1, sum2);
    shift    = norm_l(prod);
    prod     = L_shl_pos(prod, shift);
    prod_exp = sub(62, add(add(scale1, scale2), shift));
    inv      = Isqrt(prod, &prod_exp);
    scale0   = norm_l(sum0);
    sum0     = L_shl_pos(sum0, scale0);
    prod     = Mpy_32_32(sum0, inv);
    prod_exp = add(sub(31, scale0), prod_exp);
    test();
    IF (prod == 0 || sub(norm_l(prod), prod_exp) >= 0)
    {
        *normcorr = s_max(0, round_fx_sat(L_shl_sat(prod, prod_exp))); move16();
    }
    ELSE
    {
        *normcorr = 32767; move16();
    }

    /* Second try in the neighborhood of the previous pitch */
    min_pitch = s_max(MIN_PITCH_6K4, sub(*mem_pitch, 4));
    max_pitch = s_min(MAX_PITCH_6K4, add(*mem_pitch, 4));
    max32     = ac[min_pitch - MIN_PITCH_6K4]; move32();
    pitch2    = min_pitch;                     move16();
    FOR (n = min_pitch + 1; n <= max_pitch; n++)
    {
        L_tmp = L_sub_sat(ac[n - MIN_PITCH_6K4], max32);
        if (L_tmp > 0)
        {
            pitch2 = n; move16();
        }
        max32 = L_max(ac[n - MIN_PITCH_6K4], max32);
    }
    IF (sub(*pitch, pitch2) != 0)
    {
        sum0 = L_mult0(s6k4[0], s6k4[0 - pitch2]);
        sum1 = L_mac0(1, s6k4[0 - pitch2], s6k4[0 - pitch2]);
        sum2 = L_mac0(1, s6k4[0], s6k4[0]);
        for (m = 1; m < len2; m++)
        {
            sum0 = L_mac0(sum0, s6k4[m], s6k4[m - pitch2]);
            sum1 = L_mac0(sum1, s6k4[m - pitch2], s6k4[m - pitch2]);
            sum2 = L_mac0(sum2, s6k4[m], s6k4[m]);
        }
        scale1   = norm_l(sum1);
        scale2   = norm_l(sum2);
        sum1     = L_shl_pos(sum1, scale1);
        sum2     = L_shl_pos(sum2, scale2);
        prod     = Mpy_32_32(sum1, sum2);
        shift    = norm_l(prod);
        prod     = L_shl_pos(prod, shift);
        prod_exp = sub(62, add(add(scale1, scale2), shift));
        inv      = Isqrt(prod, &prod_exp);
        scale0   = norm_l(sum0);
        sum0     = L_shl_pos(sum0, scale0);
        prod     = Mpy_32_32(sum0, inv);
        prod_exp = add(sub(31, scale0), prod_exp);
        test();
        IF (prod == 0 || sub(norm_l(prod), prod_exp) >= 0)
        {
            normcorr2 = s_max(0, round_fx_sat(L_shl_sat(prod, prod_exp))); move16();
        }
        ELSE
        {
            normcorr2 = 32767; move16();
        }
        IF (sub(normcorr2, mult_r(*normcorr, 27853)) > 0)
        {
            *pitch    = pitch2;    move16();
            *normcorr = normcorr2; move16();
        }
    }
    *mem_pitch = *pitch; move16();

    /* Update memory */

    basop_memmove(mem_s6k4, &mem_s6k4[len2], MAX_PITCH_6K4 * sizeof(Word16));

    /* Upsample pitch by a factor of 2 (6.4kHz -> 12.8kHz) */
    *pitch = shl_pos(*pitch, 1); move16();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

