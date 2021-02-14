/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"

#include "constants.h"
#include "functions.h"


static Word16 getScaleFactor16_withNegativeScaling(Word16 *data16, Word16 dataLen);

void processPCapply_fx(Word16 yLen, Word16 q_old_res_fx[], Word16 *q_old_res_fx_exp, Word16 q_res_fx[],
                       Word16 q_old_d_fx[], Word16 spec_inv_idx, Word16 *fac, Word16 *fac_e,
                       Word32 q_d_fx[], Word16 *q_fx_exp, Word16 gg_idx, Word16 gg_idx_off, Word16 prev_gg, Word16 prev_gg_e,
                       Word16 *pc_nbLostFramesInRow)
{
    Counter i;
    Word16  s, s2, s3, c, tmp16, tmp16_2, inv_gain, thr;
    Word32  ener_curr, ener_prev, mean_nrg_high, mean_nrg_low;
    Word16  global_gain, global_gain_e, gg2, gg2_e, prev_gg2, prev_gg2_e;
    Word32  tmp32, ener_curr_gg2, ener_prev_gg2;

#ifdef DYNMEM_COUNT
    struct _dynmem
    {
        Counter i;
        Word16  s, s2, s3, c, tmp16, tmp16_2, inv_gain, thr;
        Word32  ener_curr, ener_prev, mean_nrg_high, mean_nrg_low;
        Word16  global_gain, global_gain_e, gg2, gg2_e, prev_gg2, prev_gg2_e;
        Word32  tmp32, ener_curr_gg2, ener_prev_gg2;
    };
    Dyn_Mem_In("processPCapply_fx", sizeof(struct _dynmem));
#endif

    assert(spec_inv_idx >= 0);

    *pc_nbLostFramesInRow = add(*pc_nbLostFramesInRow , 1);

    tmp32         = L_shl_pos(L_mult0(add(gg_idx, gg_idx_off), 0x797D), 7);
    global_gain_e = add(extract_l(L_shr_pos(tmp32, 25)), 1);
    global_gain   = round_fx(BASOP_Util_InvLog2(L_or(tmp32, 0xFE000000)));

    s = global_gain_e;  move16();
    inv_gain = Inv16(global_gain, &s);
    *fac = mult(prev_gg, inv_gain);
    *fac_e = add(s, prev_gg_e);

    /* Calculate rescaling factor */
    s = getScaleFactor16(q_old_d_fx, yLen);

    mean_nrg_low = 0;  move32();
    FOR (i = 0; i < spec_inv_idx; i++)
    {
        tmp16        = shl_sat(q_old_d_fx[i], sub(s, 4));
        mean_nrg_low = L_mac0(mean_nrg_low, tmp16, tmp16); /* exp = 2s - 8 */
    }

    mean_nrg_high = 0;  move32();
    FOR (i = spec_inv_idx; i < yLen; i++)
    {
        tmp16         = shl_sat(q_old_d_fx[i], sub(s, 4));
        mean_nrg_high = L_mac0(mean_nrg_high, tmp16, tmp16); /* exp = 2s - 8 */
    }

    IF (sub(spec_inv_idx, sub(yLen, spec_inv_idx)) < 0)
    {
        c             = div_s(spec_inv_idx, sub(yLen, spec_inv_idx));
        mean_nrg_high = Mpy_32_16(mean_nrg_high, c); /* exp = 2s - 8 */
    }
    ELSE
    {
        c            = div_s(sub(yLen, spec_inv_idx), spec_inv_idx);
        mean_nrg_low = Mpy_32_16(mean_nrg_low, c); /* exp = 2s - 8 */
    }

    IF (L_sub(mean_nrg_low, mean_nrg_high) > 0)
    {
        s         = getScaleFactor16(q_old_res_fx, spec_inv_idx);
        ener_prev = 0;  move32();
        FOR (i = 0; i < spec_inv_idx; i++)
        {
            tmp16     = shl_sat(q_old_res_fx[i], sub(s, 4));
            ener_prev = L_mac0(ener_prev, tmp16, tmp16); /* exp = - (2s - 8 - 2**q_old_res_fx_exp) */
        }

        s2        = getScaleFactor16(q_res_fx, spec_inv_idx);
        ener_curr = 0;  move32();
        FOR (i = 0; i < spec_inv_idx; i++)
        {
            tmp16     = shl_sat(q_res_fx[i], sub(s2, 4));
            ener_curr = L_mac0(ener_curr, tmp16, tmp16); /* exp = - (2s2 - 8) */
        }

        s  = shl(sub(s, *q_old_res_fx_exp), 1);
        s2 = shl(s2, 1);
        s3        = s_max(s, s2);
        ener_prev = L_shr_sat(ener_prev, sub(s3, s2));
        ener_curr = L_shr_sat(ener_curr, sub(s3,  s));

        prev_gg2   = mult(prev_gg, prev_gg);
        prev_gg2_e = shl(prev_gg_e, 1);
        ener_prev_gg2 = Mpy_32_16(ener_prev, prev_gg2); /* exp =  prev_gg2_e */

        gg2   = mult(global_gain, global_gain);
        gg2_e = shl(global_gain_e, 1);
        ener_curr_gg2 = Mpy_32_16(ener_curr, gg2);      /* exp =  gg2_e */

        s3            = s_max(prev_gg2_e, gg2_e);
        ener_prev_gg2 = L_shr_sat(ener_prev_gg2, sub(s3, prev_gg2_e));
        ener_curr_gg2 = L_shr_sat(ener_curr_gg2, sub(s3, gg2_e));

        IF (L_sub(ener_prev_gg2, ener_curr_gg2) > 0)
        {
            s       = getScaleFactor32(&ener_prev, 1);
            s2      = getScaleFactor32(&ener_curr, 1);
            s3      = s_min(s, s2);
            tmp16   = extract_h(L_shl_sat(ener_curr, s3));
            tmp16_2 = extract_h(L_shl_sat(ener_prev, s3));

            *fac_e = 0;  move16();
            tmp16_2 = Inv16(tmp16_2, fac_e);

            *fac = mult(tmp16, tmp16_2);

            IF (sub(*fac, 32767) < 0)
            {
                *fac = Sqrt16(*fac, fac_e); move16();
            }
        }
    }

    /* write synthesized samples */
    *q_old_res_fx_exp = add(*q_old_res_fx_exp, *fac_e);
    thr = shl_sat(20480, sub(-15, *q_old_res_fx_exp));
    FOR (i = spec_inv_idx; i < yLen; i++)
    {
        q_res_fx[i] = extract_h(L_mult(q_old_res_fx[i]  /* exp = q_old_res_fx_exp' */, *fac /* exp = fac_e */)); /* exp = q_old_res_fx_exp */

        IF (sub(abs_s(q_res_fx[i]), thr) < 0)
        {
            q_res_fx[i] = 0;  move16();
        }
    }

    /* scaling to 15Q16 */
    s = getScaleFactor16_withNegativeScaling(&q_res_fx[0], spec_inv_idx); /* exp = 0 */
    s2  = getScaleFactor16_withNegativeScaling(&q_res_fx[spec_inv_idx], sub(yLen, spec_inv_idx)); /* exp = q_old_res_fx_exp */
    s3  = add(s, *q_old_res_fx_exp);
    IF (sub(s3, s2) > 0) {
        tmp16 = sub(s3, s2);
        s = sub(s, tmp16);
        s3 = sub(s3, tmp16);
    }
    *q_fx_exp = sub(15, s); move16();

    FOR (i = 0; i < spec_inv_idx; i++)
    {
        q_d_fx[i] = L_shl_pos(L_deposit_h(q_res_fx[i]), s); move32();
    }
    FOR (; i < yLen; i++)
    {
        q_d_fx[i] = L_shl_pos(L_deposit_h(q_res_fx[i]), s3); move32();
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

static Word16 getScaleFactor16_withNegativeScaling(Word16 *data16, Word16 dataLen)
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

    Dyn_Mem_Deluxe_Out();

    return shift;
}


