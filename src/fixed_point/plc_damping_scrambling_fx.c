/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"
#ifdef NONBE_PLC4_ADAP_DAMP

#include "functions.h"


void processPLCDampingScrambling_main_fx(Word16 bfi, Word16 concealMethod, Word16 ns_nbLostFramesInRow,
                                         Word16 pc_nbLostFramesInRow, Word16 *ns_seed, Word16 *pc_seed, Word16 pitch_present_bfi1,
                                         Word16 pitch_present_bfi2, Word32 spec[], Word16 *q_fx_exp, Word16 *q_old_d_fx,
                                         Word16 *q_old_fx_exp, Word16 L_spec, Word16 stabFac, Word16 frame_dms,
                                         Word16 *cum_fading_slow, Word16 *cum_fading_fast, Word16 *alpha, Word16 spec_inv_idx)
{
    test();test();
    IF (( sub(bfi, 1) == 0 && sub(concealMethod, 4) == 0 ) || sub(bfi, 2) == 0)
    {
        IF (sub(bfi, 1) == 0)
        {
            processPLCDampingScrambling_fx(spec, L_spec, ns_nbLostFramesInRow, stabFac,
                                           pitch_present_bfi1, frame_dms, cum_fading_slow,
                                           cum_fading_fast, alpha, ns_seed, 0);
        }
        ELSE
        {
            processPLCDampingScrambling_fx(spec, L_spec, pc_nbLostFramesInRow, stabFac,
                                           pitch_present_bfi2, frame_dms, cum_fading_slow,
                                           cum_fading_fast, alpha, pc_seed, spec_inv_idx);

            processPLCupdateSpec_fx(q_old_d_fx, q_old_fx_exp, spec, q_fx_exp, L_spec);
        }
    }
}

void processPLCDampingScrambling_fx(Word32 spec[], Word16 L_spec, Word16 nbLostFramesInRow, Word16 stabFac,
                                    Word16 pitch_present, Word16 frame_dms, Word16 *cum_fading_slow,
                                    Word16 *cum_fading_fast, Word16 *alpha, Word16 *seed, Word16 spec_inv_idx)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
        Word16 xLostFramesInRow, slow, fast, tmp16;
        Word16 plc_start_inFrames, plc_end_inFrames, plc_duration_inFrames, x, m, b, linFuncStartStop;
        Word16 randThreshold, ad_threshFac, energThreshold, s, s2, s3, mean_energy16;
        Word32 frame_energy, mean_nrg, fac;
    );

    IF (sub(nbLostFramesInRow, 1) == 0)
    {
        *cum_fading_slow = 32767;  move16();
        *cum_fading_fast = 32767;  move16();
    }

    xLostFramesInRow = nbLostFramesInRow;  move16();
    SWITCH (frame_dms)
    {
    case 25: nbLostFramesInRow = shr(add(nbLostFramesInRow, 3), 2); BREAK;
    case 50: nbLostFramesInRow = shr(add(nbLostFramesInRow, 1), 1); BREAK;
    }

    /* get damping factor */
    tmp16 = mult(6554, stabFac);
    slow = add(26214, tmp16);
    fast = add(9830, tmp16);

    IF (spec_inv_idx == 0)
    {
        IF (sub(nbLostFramesInRow, PLC_FADEOUT_IN_MS/10) > 0)
        {
            slow = 0;  move16();
            fast = 0;  move16();
        }
        ELSE IF (sub(nbLostFramesInRow, 2) > 0)
        {
            tmp16 = div_s(sub(PLC_FADEOUT_IN_MS / 10, nbLostFramesInRow), sub(PLC_FADEOUT_IN_MS / 10, sub(nbLostFramesInRow, 1)));
            slow = mult(slow, tmp16);
            fast = mult(fast, tmp16);
        }
    }

    SWITCH (frame_dms)
    {
    case 25:
        IF (sub(slow, 32767) < 0)
        {
            tmp16  = 0;
            slow = Sqrt16(slow, &tmp16);  move16();
            slow = shl(slow, tmp16);
        }
        IF (sub(slow, 32767) < 0)
        {
            tmp16  = 0;
            slow = Sqrt16(slow, &tmp16);  move16();
            slow = shl(slow, tmp16);
        }
        IF (sub(fast, 32767) < 0)
        {
            tmp16  = 0;
            fast = Sqrt16(fast, &tmp16);  move16();
            fast = shl(fast, tmp16);
        }
        IF (sub(fast, 32767) < 0)
        {
            tmp16  = 0;
            fast = Sqrt16(fast, &tmp16);  move16();
            fast = shl(fast, tmp16);
        }
        BREAK;
    case 50:
        IF (sub(slow, 32767) < 0)
        {
            tmp16  = 0;
            slow = Sqrt16(slow, &tmp16);  move16();
            slow = shl(slow, tmp16);
        }
        IF (sub(fast, 32767) < 0)
        {
            tmp16  = 0;
            fast = Sqrt16(fast, &tmp16);  move16();
            fast = shl(fast, tmp16);
        }
        BREAK;
    }

    *alpha = slow;
    *cum_fading_slow = mult_r(*cum_fading_slow, slow);
    *cum_fading_fast = mult_r(*cum_fading_fast, fast);

    /* Get fadeout function */
    /* being 1 up to plc_start_inFrames, being 0 starting with
       plc_end_inFrames; decreasing linearly in between */
    SWITCH (frame_dms)
    {
    case 25:
        plc_start_inFrames = (10*PLC_START_IN_MS) / 25;  move16();
        plc_end_inFrames   = (10*PLC_FADEOUT_IN_MS) / 25;  move16();
        BREAK;
    case 50:
        plc_start_inFrames = PLC_START_IN_MS / 5;  move16();
        plc_end_inFrames   = PLC_FADEOUT_IN_MS / 5;  move16();
        BREAK;
    default:
        plc_start_inFrames = PLC_START_IN_MS / 10;  move16();
        plc_end_inFrames   = PLC_FADEOUT_IN_MS / 10;  move16();
    }

    IF (pitch_present == 0)
    {
        plc_start_inFrames = 1;  move16();
    }
    plc_duration_inFrames = sub(plc_end_inFrames, plc_start_inFrames);

    IF (sub(xLostFramesInRow, plc_start_inFrames) < 0)
    {
        x = plc_start_inFrames;  move16();
    }
    ELSE
    {
        if (sub(xLostFramesInRow, plc_end_inFrames) > 0)
        {
            x = plc_end_inFrames;  move16();
        }
        else
        {
            x = xLostFramesInRow;  move16();
        }
    }
    m = negate(div_s(1, plc_duration_inFrames));
    b = negate(plc_end_inFrames); /* shift on x axis */
    linFuncStartStop = i_mult(add(m, 1), (add(x, b)));  /* +1 to prevent 32768 as result */

    /* Sign scrambling */
    randThreshold = mult(-32768, linFuncStartStop);

    tmp16 = *seed;  move16();
    FOR (i = spec_inv_idx; i < L_spec; i++)
    {
        tmp16 = extract_l(L_mac0(16831, tmp16, 12821));

        IF (tmp16 < 0)
        {
            test();
            if (pitch_present == 0 || sub(tmp16, randThreshold) < 0)
            {
                spec[i] = L_negate(spec[i]);
            }
        }

    }
    *seed = tmp16; move16();

    /* Apply adaptive damping */
    tmp16 = mult(18022 /* 10 - 1.2 */, linFuncStartStop);
    ad_threshFac = add(shr(tmp16, 1), 1228 /* 1.2 >> 1 */); /* exp = 5 */

    s = getScaleFactor32(&spec[spec_inv_idx], sub(L_spec, spec_inv_idx));
    frame_energy = 0;  move32();
    FOR (i = spec_inv_idx; i < L_spec; i++)
    {
        tmp16     = extract_h(L_shl_sat(spec[i], sub(s, 4)));
        frame_energy = L_mac0(frame_energy, tmp16, tmp16); /* exp = -(2*(s-16) - 8) */
    }
    mean_energy16 = BASOP_Util_Divide3216_Scale(frame_energy, sub(L_spec, spec_inv_idx), &s2);  /* exp = -(2*(s-16) - 8) + 16 - (15-s2) */

    energThreshold = mult(ad_threshFac, mean_energy16);    /* exp = -(2*(s-16) - 8) + 16 - (15-s2) + 5 */

    s3 = add(sub(29, shl(sub(s, 16), 1)), s2);
    IF (sub(energThreshold, 32767) < 0)
    {
        energThreshold = Sqrt16(energThreshold, &s3);
    }
    s3 = sub(s3, 15);

    mean_nrg = L_shl_sat(L_deposit_l(energThreshold), s3); /* exp = 0 */
    fac = mult(sub(*cum_fading_slow, *cum_fading_fast), energThreshold);
    fac = L_shl_sat(L_deposit_l(fac), s3); /* exp = 0 */

    FOR (i = spec_inv_idx; i < L_spec; i++)
    {
        if (L_sub(L_abs(spec[i]), mean_nrg) < 0)
        {
            spec[i] = Mpy_32_16(spec[i], *cum_fading_slow);
        }
        else
        {
            if (spec[i] > 0)
            {
                spec[i] = L_add(Mpy_32_16(spec[i], *cum_fading_fast), fac);
            }
            else if (spec[i] == 0)
            {
                spec[i] = Mpy_32_16(spec[i], *cum_fading_fast);
            }
            else
            {
                spec[i] = L_sub(Mpy_32_16(spec[i], *cum_fading_fast), fac);
            }
        }
    }

    Dyn_Mem_Deluxe_Out();
}


#endif /* NONBE_PLC4_ADAP_DAMP */
