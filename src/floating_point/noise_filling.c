/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processNoiseFilling_fl(LC3_FLOAT xq[], LC3_INT nfseed, LC3_INT fac_ns_idx, LC3_INT BW_cutoff_idx, LC3_INT frame_dms)
{
    LC3_INT   zeroLines[MAX_LEN] = {0};
    LC3_INT   nTransWidth = 0, startOffset = 0, i = 0, j = 0, k = 0, start = 0, end = 0, allZeros = 0, kZeroLines = 0;
    LC3_FLOAT fac_ns = 0;

    switch (frame_dms)
    {
        case 25:
            nTransWidth = 1;
            startOffset = 6;
            break;
        case 50:
            nTransWidth = 1;
            startOffset = 12;
            break;
        case 100:
            nTransWidth = 3;
            startOffset = 24;
            break;
    }

    fac_ns = (8.0 - fac_ns_idx) / 16.0;

    j = 0;

    for (k = startOffset; k < BW_cutoff_idx; k++) {
        allZeros = 1;

        start = k - nTransWidth;
        end   = MIN(BW_cutoff_idx - 1, k + nTransWidth);

        for (i = start; i <= end; i++) {
            if (xq[i] != 0) {
                allZeros = 0;
            }
        }

        if (allZeros == 1) {
            zeroLines[j] = k;
            kZeroLines++;
            j++;
        }
    }

    for (k = 0; k < kZeroLines; k++) {
        nfseed = (13849 + (nfseed + 32768) * 31821) & 65535;
        nfseed -= 32768;

        if (nfseed >= 0) {
            xq[zeroLines[k]] = fac_ns;
        } else {
            xq[zeroLines[k]] = -fac_ns;
        }
    }

    return;
}
