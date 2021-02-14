/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"



void processAdjustGlobalGain_fx(Word16 *gg_idx, Word16 gg_idx_min, Word16 gg_idx_off, Word16 *gain, Word16 *gain_e,
                                Word16 target, Word16 nBits, Word16 *gainChange, Word16 fs_idx)
{

    Word32 L_tmp;
    Word16 delta, delta2;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processAdjustGlobalGain_fx", sizeof(struct {
                   Word32 L_tmp;
                   Word16 delta, delta2;
               }));
#endif

    IF (sub(nBits, adjust_global_gain_tables[0][fs_idx]) < 0)
    {
        delta = mult_r(add(nBits, 48), 2048);
    }
    ELSE IF (sub(nBits, adjust_global_gain_tables[1][fs_idx]) < 0)
    {
        delta = mult_r(add(nBits, adjust_global_gain_tables[4][fs_idx]), adjust_global_gain_tables[3][fs_idx]);
    }
    ELSE IF (sub(nBits, adjust_global_gain_tables[2][fs_idx]) < 0)
    {
        delta = mult_r(nBits, 683);
    }
    ELSE
    {
        delta = mult_r(adjust_global_gain_tables[2][fs_idx], 683);
    }
    delta2 = add(delta, 2);

    *gainChange = 0; move16();

    test();
    IF (sub(*gg_idx, 255) == 0 && sub(nBits, target) > 0)
    {
        *gainChange = 1; move16();
    }

    test(); test(); test();
    IF ((sub(*gg_idx, 255) < 0 && sub(nBits, target) > 0) || (*gg_idx > 0 && sub(nBits, sub(target, delta2)) < 0))
    {
        test();
        IF (sub(nBits, sub(target, delta2)) < 0)
        {
            *gg_idx = sub(*gg_idx, 1); move16();
        }
        ELSE IF (sub(*gg_idx, 254) == 0 || sub(nBits, add(target, delta)) < 0)
        {
            *gg_idx = add(*gg_idx, 1); move16();
        }
        ELSE
        {
            *gg_idx = add(*gg_idx, 2); move16();
        }

        *gg_idx = s_max(*gg_idx, sub(gg_idx_min, gg_idx_off)); move16();

        L_tmp       = L_shl_pos(L_mult0(add(*gg_idx, gg_idx_off), 0x797D), 7); /* 6Q25; 0x797D -> log2(10)/28 (Q18) */
        *gain_e     = add(extract_l(L_shr_pos(L_tmp, 25)), 1);                 /* get exponent */
        *gain       = round_fx(BASOP_Util_InvLog2(L_or(L_tmp, 0xFE000000)));
        *gainChange = 1; move16();
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

