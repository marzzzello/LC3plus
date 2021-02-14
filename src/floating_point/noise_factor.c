/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processNoiseFactor_fl(LC3_INT* fac_ns_idx, LC3_FLOAT x[], LC3_INT xq[], LC3_FLOAT gg, LC3_INT BW_cutoff_idx, LC3_INT frame_dms,
                           LC3_INT target_bytes
)
{
    LC3_INT sumZeroLines = 0, kZeroLines = 0, startOffset = 0, nTransWidth = 0, end = 0, start = 0, i = 0, j = 0, k = 0,
        allZeros = 0, m = 0;
    LC3_FLOAT fac_ns_unq = 0, mean = 0, idx = 0, nsf1 = 0, nsf2 = 0;
    LC3_INT   zeroLines[MAX_LEN] = {0}, zL1[MAX_LEN] = {0}, zL2[MAX_LEN] = {0};

    switch (frame_dms)
    {
        case 25:
            nTransWidth = 4;
            startOffset = 6;
            break;
        case 50:
            nTransWidth = 4;
            startOffset = 12;
            break;
        case 100:
            nTransWidth = 8;
            startOffset = 24;
            break;
    }

    for (k = startOffset; k < BW_cutoff_idx; k++) {
        allZeros = 1;

        start = k - (nTransWidth - 2) / 2;
        end   = MIN(BW_cutoff_idx - 1, k + (nTransWidth - 2) / 2);

        for (i = start; i <= end; i++) {
            if (xq[i] != 0) {
                allZeros = 0;
            }
        }

        if (allZeros == 1) {
            zeroLines[j] = k + 1;
            kZeroLines++;
            j++;
        }
    }

    for (i = 0; i < kZeroLines; i++) {
        sumZeroLines += zeroLines[i];
    }

    if (sumZeroLines > 0) {
        for (j = 0; j < kZeroLines; j++) {
            mean += LC3_FABS(x[zeroLines[j] - 1] / gg);
        }

        fac_ns_unq = mean / kZeroLines;
    } else {
        fac_ns_unq = 0;
    }

    if (target_bytes <= 20 && frame_dms == 100) {
        j = 0, k = 0;
        m = floor(sumZeroLines / kZeroLines);

        for (i = 0; i < kZeroLines; i++) {
            if (zeroLines[i] <= m) {
                zL1[j] = zeroLines[i];
                j++;
            }

            if (zeroLines[i] > m) {
                zL2[k] = zeroLines[i];
                k++;
            }
        }

        mean = 0;
        for (i = 0; i < j; i++) {
            mean += LC3_FABS(x[zL1[i] - 1]) / gg;
        }

        nsf1 = mean / j;

        mean = 0;
        for (i = 0; i < k; i++) {
            mean += LC3_FABS(x[zL2[i] - 1]) / gg;
        }

        nsf2 = mean / k;

        fac_ns_unq = MIN(nsf1, nsf2);
    }

    idx = round(8 - 16 * fac_ns_unq);
    idx = MIN(MAX(idx, 0), 7);

    *fac_ns_idx = idx;
}
