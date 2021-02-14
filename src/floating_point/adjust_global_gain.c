/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processAdjustGlobalGain_fl(LC3_INT* gg_idx, LC3_INT gg_idx_min, LC3_INT gg_idx_off, LC3_FLOAT* gain, LC3_INT target, LC3_INT nBits, LC3_INT* gainChange, LC3_INT fs_idx)
{
    LC3_FLOAT delta  = 0;
    LC3_INT   delta2 = 0;

    if (nBits < gg_p1[fs_idx]) {
        delta = (nBits + 48.0) / 16.0;
    } else if (nBits < gg_p2[fs_idx]) {
        delta = (nBits + gg_d[fs_idx]) * gg_c[fs_idx];
    } else if (nBits < gg_p3[fs_idx]) {
        delta = nBits / 48.0;
    } else {
        delta = gg_p3[fs_idx] / 48.0;
    }

    delta  = round(delta);
    delta2 = delta + 2;

    *gainChange = 0;

    if (*gg_idx == 255 && nBits > target) {
        *gainChange = 1;
    }

    if ((*gg_idx < 255 && nBits > target) || (*gg_idx > 0 && nBits < target - delta2)) {
        if (nBits < target - delta2) {
            *gg_idx = *gg_idx - 1;
        } else if (*gg_idx == 254 || nBits < target + delta) {
            *gg_idx = *gg_idx + 1;
        } else {
            *gg_idx = *gg_idx + 2;
        }

        *gg_idx     = MAX(*gg_idx, gg_idx_min - gg_idx_off);
        *gain       = LC3_POW(10, (LC3_FLOAT)(*gg_idx + gg_idx_off) / 28);
        *gainChange = 1;
    }
}
