/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processTnsDecoder_fl(LC3_FLOAT* x, LC3_INT* rc_idx, LC3_INT* order, LC3_INT numfilters, LC3_INT bw_fcbin, LC3_INT N, LC3_INT fs, LC3_INT frame_dms)
{
    LC3_INT   startfreq[2] = {0}, stopfreq[2] = {0}, f = 0, i = 0, j = 0, m = 0, l = 0, rc_idx_f[9] = {0};
    LC3_FLOAT rc[9] = {0}, s = 0, st[9] = {0};

    if (N > 40 * ((LC3_FLOAT) (frame_dms) / 10.0)) {
        N  = 40 * ((LC3_FLOAT) (frame_dms) / 10.0);
        fs = 40000;
    }

    if (numfilters == 2) {
        startfreq[0] = floor(600 * N * 2 / fs) + 1;
        stopfreq[0]  = bw_fcbin / 2;
        startfreq[1] = bw_fcbin / 2 + 1;
        stopfreq[1]  = bw_fcbin;
    } else {
        startfreq[0] = floor(600 * N * 2 / fs) + 1;
        stopfreq[0]  = bw_fcbin;
    }

    for (f = 0; f < numfilters; f++) {
        if (order[f] > 0) {
            j = 0;

            for (i = f * 8; i < f * 8 + 8; i++) {
                rc_idx_f[j] = rc_idx[i];
                rc[j]       = quants_pts_tns[rc_idx_f[j]];
                j++;
            }

            for (m = startfreq[f]; m <= stopfreq[f]; m++) {
                s = x[m - 1] - rc[order[f] - 1] * st[order[f] - 1];

                for (l = order[f] - 2; l >= 0; l--) {
                    s         = s - rc[l] * st[l];
                    st[l + 1] = rc[l] * s + st[l];
                }

                st[0]    = s;
                x[m - 1] = s;
            }
        }
    }
}
