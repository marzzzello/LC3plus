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


void processPCupdate_fx(Word16 bfi, Word16 yLen, Word16 q_old_res_fx[], Word16 *q_old_res_fx_exp,
                      Word16 q_res_fx[], Word16 spec_inv_idx, Word16 gg_idx, Word16 gg_idx_off,
                      Word16 *prev_gg, Word16 *prev_gg_e, Word16 rframe, Word16 *BW_cutoff_idx_nf,
                      Word16 *prev_BW_cutoff_idx_nf, Word16 fac_ns_idx, Word16 *prev_fac_ns_fx, Word16 fac,
                      Word16 fac_e)
{
    Word16  global_gain, global_gain_e, s, s2, s3, tmp16;
    Word32  tmp32;

#ifdef DYNMEM_COUNT
    struct _dynmem
    {
        Word16  global_gain, global_gain_e, s, s2, s3, tmp16;
        Word32  tmp32;
    };
    Dyn_Mem_In("processPCupdate_fx", sizeof(struct _dynmem));
#endif

    tmp32         = L_shl_pos(L_mult0(add(gg_idx, gg_idx_off), 0x797D), 7);
    global_gain_e = add(extract_l(L_shr_pos(tmp32, 25)), 1);
    global_gain = round_fx(BASOP_Util_InvLog2(L_or(tmp32, 0xFE000000)));

    *prev_gg   = global_gain;  move16();
    *prev_gg_e = global_gain_e;  move16();

    s = getScaleFactor16(q_res_fx, spec_inv_idx); /* exp = 0 */
    IF (bfi == 0)
    {
        *q_old_res_fx_exp = negate(s);
        Copy_Scale_sig(q_res_fx, q_old_res_fx, yLen, s);
    }
    ELSE
    {
        s2 = getScaleFactor16(&q_res_fx[spec_inv_idx], sub(yLen, spec_inv_idx)); /* exp = q_old_res_fx_exp */
        s3 = add(s, *q_old_res_fx_exp);
        IF (sub(s3, s2) > 0)
        {
            tmp16 = sub(s3, s2);
            s = sub(s, tmp16);
        }
        s2 = add(s, *q_old_res_fx_exp);
        *q_old_res_fx_exp = negate(s);

        Copy_Scale_sig(q_res_fx, q_old_res_fx, spec_inv_idx, s);
        Copy_Scale_sig(&q_res_fx[spec_inv_idx], &q_old_res_fx[spec_inv_idx], sub(yLen, spec_inv_idx), s2);
    }

    IF (rframe == 0)
    {
        *prev_BW_cutoff_idx_nf = *BW_cutoff_idx_nf;
        *prev_fac_ns_fx = shl_pos(sub(8, fac_ns_idx), 11);
    }
    ELSE IF(sub(bfi, 2) == 0 && sub(*BW_cutoff_idx_nf, *prev_BW_cutoff_idx_nf) != 0)
    {
        *BW_cutoff_idx_nf = *prev_BW_cutoff_idx_nf;
        *prev_fac_ns_fx = shl_sat(mult(*prev_fac_ns_fx, fac), fac_e);
        *prev_fac_ns_fx = s_max(*prev_fac_ns_fx, 2048);
        *prev_fac_ns_fx = s_min(*prev_fac_ns_fx, 16384);
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}


