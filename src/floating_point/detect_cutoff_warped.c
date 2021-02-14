/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processDetectCutoffWarped_fl(LC3_FLOAT* d2, LC3_INT fs_idx, LC3_INT frame_dms, LC3_INT* bw_idx)
{
    const LC3_INT *warp_idx_start = NULL, *warp_idx_stop = NULL;
    LC3_INT        counter = 0, brickwall = 0, i = 0, stop = 0, dist = 0;
    LC3_FLOAT      d2_mean = 0, d2_sum = 0, e_diff = 0, thr = 0;

    warp_idx_start = BW_warp_idx_start_all[fs_idx - 1];
    warp_idx_stop  = BW_warp_idx_stop_all[fs_idx - 1];

    switch (frame_dms)
    {
        case 25:
            warp_idx_start = BW_warp_idx_start_all_2_5ms[fs_idx - 1];
            warp_idx_stop  = BW_warp_idx_stop_all_2_5ms[fs_idx - 1];
            break;
        case 50:
            warp_idx_start = BW_warp_idx_start_all_5ms[fs_idx - 1];
            warp_idx_stop  = BW_warp_idx_stop_all_5ms[fs_idx - 1];
            break;
        case 100:
            warp_idx_start = BW_warp_idx_start_all[fs_idx - 1];
            warp_idx_stop  = BW_warp_idx_stop_all[fs_idx - 1];
            break;
    }
    
    counter = fs_idx;
    
    d2_sum = sum_vec(&d2[warp_idx_start[counter - 1]], warp_idx_stop[counter - 1] - warp_idx_start[counter - 1] + 1);

    d2_mean = d2_sum / (warp_idx_stop[counter - 1] - warp_idx_start[counter - 1] + 1);

    while (d2_mean < threshold_quiet[counter - 1]) {
        d2_sum = 0;
        counter--;
        if (counter == 0) {
            break;
        }

        /* calculate mean energy per band */
        d2_sum =
            sum_vec(&d2[warp_idx_start[counter - 1]], warp_idx_stop[counter - 1] - warp_idx_start[counter - 1] + 1);

        d2_mean = d2_sum / (warp_idx_stop[counter - 1] - warp_idx_start[counter - 1] + 1);
    }

    *bw_idx = counter;

    /* check if energy difference between bands is present */
    if (*bw_idx < fs_idx) {
        thr  = (LC3_FLOAT)threshold_brickwall[counter];
        stop = warp_idx_start[counter];
        dist = brickwall_dist[counter];

        for (i = stop; i >= stop - dist; i--) {
            e_diff = 10.0 * LC3_LOG10(d2[i - dist + 1] + FLT_EPSILON) - 10.0 * LC3_LOG10(d2[i + 1] + FLT_EPSILON);

            if (e_diff > thr) {
                brickwall = 1;
                break;
            }
        }

        if (brickwall == 0) {
            *bw_idx = fs_idx;
        }
    }
}
