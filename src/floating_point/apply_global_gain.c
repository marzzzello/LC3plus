/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processApplyGlobalGain_fl(LC3_FLOAT x[], LC3_INT xLen, LC3_INT global_gain_idx, LC3_INT global_gain_off)
{
    LC3_FLOAT gg = 0;

    gg = LC3_POW(10, (LC3_FLOAT)(global_gain_idx + global_gain_off) / 28);

    mult_const(x, gg, xLen);
}
