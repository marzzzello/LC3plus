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


void process_ltpf_coder_fx(Word16 *bits, Word16 ol_pitch, Word16 ltpf_enable, Word16 *mem_in_exp, Word16 mem_in[],
                           Word16 mem_in_len, Word16 param[], Word16 *xin, Word16 len, Word16 *mem_normcorr,
                           Word16 *mem_mem_normcorr, Word16 ol_normcorr, Word16 *mem_ltpf_on, Word16 *mem_ltpf_pitch,
                           Word16 xin_exp, Word16 frame_dms, Word8 *scratchBuffer)
{
    Dyn_Mem_Deluxe_In(
        Word16  pitch_index, scale0, scale1, scale2, s, *x, x_exp, shift, prod_exp, ltpf_pitch;
        Word32 *ac32;
        Word16 *ac, *currFrame, *predFrame;
        Word32  L_tmp, ac32_max, cor_max32, sum0, sum1, sum2, prod, inv;
        Word16  min_pitch, max_pitch, ac_min_pitch, ac_max_pitch, ac_max;
        Word16  pitch, pitch_res, min_pitch_fr, pitch_int, pitch_fr, norm_corr, ltpf_active;
        Word32  sum;
        Counter n, m, fr;
        Word16  tmp;
    );



    ac32      = (Word32 *)scratchAlign(scratchBuffer, 0);                         /* Size = 4 * 17 = 68 bytes;   */
    ac        = (Word16 *)scratchAlign(ac32, sizeof(*ac32) * 17);                 /* Size = 2 * 17 = 34 bytes    */
    currFrame = (Word16 *)scratchAlign(scratchBuffer, 0);                         /* Size = 2 * 128 = 256 bytes  */
    predFrame = (Word16 *)scratchAlign(currFrame, sizeof(*currFrame) * LEN_12K8); /* Size = 2 * 128 = 256 bytes  */
    /* Buffers 'overlap' since they are not used at the same time */              /* Total size used = 512 bytes */

    ltpf_active = 0; move16();
    norm_corr   = 0; move16();

    /* Input buffer */
    x = mem_in + mem_in_len;

    basop_memmove(x, xin, (len + 1) * sizeof(Word16));

    ASSERT(mem_in_len + len + 1 <= LTPF_MEMIN_LEN + LEN_12K8 + 1);

    /* Scaling */
    scale0      = sub(getScaleFactor16_0(mem_in, mem_in_len), 3);
    *mem_in_exp = sub(*mem_in_exp, scale0); move16();
    scale1      = sub(getScaleFactor16_0(x, len + 1), 3);
    x_exp       = sub(xin_exp, scale1);
    scale2      = sub(*mem_in_exp, x_exp);
    IF (scale2 > 0)
    {
        Scale_sig(x, len + 1, sub(scale1, scale2));
        Scale_sig(mem_in, mem_in_len, scale0);
        x_exp = *mem_in_exp; move16();
    }
    ELSE
    {
        Scale_sig(x, len + 1, scale1);
        Scale_sig(mem_in, mem_in_len, add(scale0, scale2));
        *mem_in_exp = x_exp; move16();
    }

    IF (sub(ol_normcorr, 19660) > 0)
    {
        /* Autocorrelation Bounds */
        min_pitch    = sub(ol_pitch, 4);
        max_pitch    = add(ol_pitch, 4);
        min_pitch    = s_max(min_pitch, MIN_PITCH_12K8);
        max_pitch    = s_min(max_pitch, MAX_PITCH_12K8);
        ac_min_pitch = sub(min_pitch, 4);
        ac_max_pitch = add(max_pitch, 4);

        /* Compute Autocorrelation */
        ac32_max = L_deposit_l(0);
        FOR (n = ac_min_pitch; n <= ac_max_pitch; n++)
        {
            sum = L_mac0(0L, x[0], x[0 - n]);
            FOR (m = 1; m < len; m++)
            {
                sum = L_mac0(sum, x[m], x[m - n]);
            }
            ac32[n - ac_min_pitch] = sum; move32();
            ac32_max               = L_max(L_abs(ac32[n - ac_min_pitch]), ac32_max);
        }
        s = norm_l(ac32_max);
        FOR (n = ac_min_pitch; n <= ac_max_pitch; n++)
        {
            ac[n - ac_min_pitch] = round_fx_sat(L_shl_sat(ac32[n - ac_min_pitch], s)); move16();
        }

        /* Find maximum */
        ac_max = ac[min_pitch - ac_min_pitch]; move16();
        pitch  = min_pitch;                    move16();
        FOR (n = min_pitch + 1; n <= max_pitch; n++)
        {
            tmp = sub_sat(ac[n - ac_min_pitch], ac_max);
            if (tmp > 0)
            {
                pitch = n; move16();
            }
            ac_max = s_max(ac_max, ac[n - ac_min_pitch]);
        }
        pitch_int   = pitch; move16();
        pitch_fr    = 0;     move16();
        pitch_index = add(pitch_int, 283);

        /* If the pitch is low -> estimate a fractional part */
        IF (sub(pitch, RES2_PITCH_12K8) < 0)
        {
            IF (sub(pitch, RES4_PITCH_12K8) < 0)
            {
                pitch_res    = 1;  move16();
                min_pitch_fr = -3; move16();
            }
            ELSE
            {
                pitch_res    = 2;  move16();
                min_pitch_fr = -2; move16();
            }
            if (sub(pitch, min_pitch) == 0)
            {
                min_pitch_fr = 0;
            }
            cor_max32 = MIN_32;
            FOR (fr = min_pitch_fr; fr < 4; fr += pitch_res)
            {
                sum = L_mult0(ac[pitch_int - ac_min_pitch - 4], ltpf_ac_interp_filt[fr + 3][0]);
                sum = L_mac0(sum, ac[pitch_int - ac_min_pitch - 3], ltpf_ac_interp_filt[fr + 3][1]);
                sum = L_mac0(sum, ac[pitch_int - ac_min_pitch - 2], ltpf_ac_interp_filt[fr + 3][2]);
                sum = L_mac0(sum, ac[pitch_int - ac_min_pitch - 1], ltpf_ac_interp_filt[fr + 3][3]);
                sum = L_mac0(sum, ac[pitch_int - ac_min_pitch + 0], ltpf_ac_interp_filt[fr + 3][4]);
                sum = L_mac0(sum, ac[pitch_int - ac_min_pitch + 1], ltpf_ac_interp_filt[fr + 3][5]);
                sum = L_mac0(sum, ac[pitch_int - ac_min_pitch + 2], ltpf_ac_interp_filt[fr + 3][6]);
                sum = L_mac0(sum, ac[pitch_int - ac_min_pitch + 3], ltpf_ac_interp_filt[fr + 3][7]);
                sum = L_mac0(sum, ac[pitch_int - ac_min_pitch + 4], ltpf_ac_interp_filt[fr + 3][8]);

                L_tmp = L_sub_sat(sum, cor_max32);
                if (L_tmp > 0)
                {
                    pitch_fr = fr; move16();
                }
                cor_max32 = L_max(cor_max32, sum);
            }
            IF (pitch_fr < 0)
            {
                pitch_int = sub(pitch_int, 1);
                pitch_fr  = add(pitch_fr, 4);
            }
            IF (sub(pitch_int, 127) >= 0)
            {
                pitch_index = add(add(shl_pos(pitch_int, 1), shr_pos(pitch_fr, 1)), 126);
            }
            ELSE
            {
                pitch_index = sub(add(shl_pos(pitch_int, 2), pitch_fr), 128);
            }
        }
        ltpf_pitch = add(shl_pos(pitch_int, 2), pitch_fr);

        /* Filter current and predicted frame */
#ifndef NONBE_PLC_PITCH_TUNING
        IF (sub(ltpf_enable, 1) == 0)
        {
#endif
        SWITCH (frame_dms)
        {
        case 25: x -= add(len, add(len, len)); BREAK; /* (3 * len) */
        case 50: x -= len; BREAK;
        case 100: BREAK;
        }

        FOR (n = 0; n < LEN_12K8; n++)
        {
            sum          = L_mult(x[n + 1], inter_filter[0][0][0]);
            sum          = L_mac(sum, x[n], inter_filter[0][0][1]);
            currFrame[n] = mac_r(sum, x[n - 1], inter_filter[0][0][2]);

            sum          = L_mult(x[n - pitch_int + 1], inter_filter[0][pitch_fr][0]);
            sum          = L_mac(sum, x[n - pitch_int], inter_filter[0][pitch_fr][1]);
            sum          = L_mac(sum, x[n - pitch_int - 1], inter_filter[0][pitch_fr][2]);
            predFrame[n] = mac_r(sum, x[n - pitch_int - 2], inter_filter[0][pitch_fr][3]);
        }

        /* Normalized Correlation */
        sum0 = L_mult0(currFrame[0], predFrame[0]);
        sum1 = L_mac0(1, predFrame[0], predFrame[0]);
        sum2 = L_mac0(1, currFrame[0], currFrame[0]);
        for (m = 1; m < LEN_12K8; m++)
        {
            sum0 = L_mac0(sum0, currFrame[m], predFrame[m]);
            sum1 = L_mac0(sum1, predFrame[m], predFrame[m]);
            sum2 = L_mac0(sum2, currFrame[m], currFrame[m]);
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
            norm_corr = s_max(0, round_fx_sat(L_shl_sat(prod, prod_exp))); move16();
        }
        ELSE
        {
            norm_corr = 32767; move16();
        }
        if (norm_corr < 0)
        {
            norm_corr = 0;
        }

#ifdef NONBE_PLC_PITCH_TUNING
        IF (sub(ltpf_enable, 1) == 0)
        {
#endif
            test(); test(); test(); test();
            /* Decision if lptf active */
            IF ((*mem_ltpf_on == 0 && sub(*mem_normcorr, 30802) > 0 && sub(norm_corr, 30802) > 0 &&
                 (sub(frame_dms, 100) == 0 || sub(*mem_mem_normcorr, 30802) > 0)) ||
                (sub(*mem_ltpf_on, 1) == 0 && sub(norm_corr, 29491) > 0) ||
                (sub(*mem_ltpf_on, 1) == 0 && sub(abs_s(sub(ltpf_pitch, *mem_ltpf_pitch)), 8) < 0 &&
                 add(sub(norm_corr, *mem_normcorr), 3277) > 0 && sub(norm_corr, 27525) > 0))
            {
                ltpf_active = 1; move16();
            }
        }

#ifdef NONBE_PLC_PITCH_TUNING
        test();
        IF (sub(frame_dms, 100) < 0 && sub(norm_corr, 22282) < 0)
        {
            param[0] = 0;  move16();
            param[1] = 0;  move16();
            param[2] = 0;  move16();
            *bits    = 1;  move16();
        }
        ELSE
        {
#endif
            param[0] = 1;           move16();
            param[1] = ltpf_active; move16();
            param[2] = pitch_index; move16();
            *bits    = 11;          move16();
#ifdef NONBE_PLC_PITCH_TUNING
        }
#endif
    }
    ELSE
    {
        norm_corr  = ol_normcorr; move16();
        param[0]   = 0;           move16();
        param[1]   = 0;           move16();
        param[2]   = 0;           move16();
        *bits      = 1;           move16();
        ltpf_pitch = 0;           move16();
    }

/* Update memory */
    FOR (n = 0; n < mem_in_len; n++)
    {
        mem_in[n] = mem_in[n + len]; move16();
    }

    *mem_mem_normcorr = *mem_normcorr; move16();
    *mem_normcorr   = norm_corr;   move16();
    *mem_ltpf_on    = ltpf_active; move16();
    *mem_ltpf_pitch = ltpf_pitch;  move16();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

