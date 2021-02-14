/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "functions.h"

void processDetectCutoffWarped_fx(Word16 *bw_idx, Word32 *d2_fx, Word16 d2_fx_exp, Word16 fs_idx, Word16 frame_dms)
{

    Dyn_Mem_Deluxe_In(
        Counter       iBand;
        Word32        d2_fx_sum;
        Word32        d2_fx_mean;
        Word32        delta_energy;
        Word16        d2_fx_sum_exp;
        Word16        d2_fx_mean_exp;
        Word16        nrg_below_thresh;
        Word16        counter;
        Word16        brickwall;
        Word16        stop;
        Word16        brickwall_dist;
        const Word16 *warp_idx_start, *warp_idx_stop;
    );

    SWITCH (frame_dms)
    {
    case 25:
        warp_idx_start = BW_warp_idx_start_all_2_5ms[fs_idx - 1]; move16();
        warp_idx_stop  = BW_warp_idx_stop_all_2_5ms[fs_idx - 1];  move16();
        BREAK;
    case 50:
        warp_idx_start = BW_warp_idx_start_all_5ms[fs_idx - 1]; move16();
        warp_idx_stop  = BW_warp_idx_stop_all_5ms[fs_idx - 1];  move16();
        BREAK;
    default:                                                /* 100 */
        warp_idx_start = BW_warp_idx_start_all[fs_idx - 1]; move16();
        warp_idx_stop  = BW_warp_idx_stop_all[fs_idx - 1];  move16();
        BREAK;
    }

    counter = fs_idx;
    DO
    {

        /* counter is 0...num_idxs-1 */
        counter = sub(counter, 1);

        /* always code the lowest band (NB), skip check against threshold if counter == -1 */
        IF (counter < 0)
        {
            BREAK;
        }

        d2_fx_mean     = 0; move32();
        d2_fx_mean_exp = 0; move16();

        iBand         = warp_idx_start[counter]; move16();
        d2_fx_sum     = d2_fx[iBand];            move32();
        d2_fx_sum_exp = d2_fx_exp;               move16();

        iBand++;
        FOR (; iBand <= warp_idx_stop[counter]; iBand++)
        {
            d2_fx_sum = BASOP_Util_Add_Mant32Exp(d2_fx[iBand], d2_fx_exp, d2_fx_sum, d2_fx_sum_exp, &d2_fx_sum_exp);
        }
        /* Energy-sum */
        d2_fx_mean = Mpy_32_16(d2_fx_sum, InvIntTable[add(sub(warp_idx_stop[counter], warp_idx_start[counter]), 1)]);
        d2_fx_mean_exp = d2_fx_sum_exp; move16();

        /* check if above threshold */
        nrg_below_thresh = BASOP_Util_Cmp_Mant32Exp(BW_thresh_quiet[counter], BW_thresh_quiet_exp, d2_fx_mean,
                                                    d2_fx_mean_exp); /* true if firstNumber > secondNumber */
    }
    WHILE (nrg_below_thresh > 0)
        ;

    *bw_idx = add(1, counter); move16();

    /* addtional check for brickwall characteristic */
    IF (sub(fs_idx, *bw_idx) > 0)
    {
        brickwall      = 0; move16();
        stop           = add(warp_idx_start[counter + 1], 1);
        brickwall_dist = BW_brickwall_dist[counter + 1];

        FOR (iBand = stop; iBand >= sub(stop, brickwall_dist); iBand--)
        {
            /* Band(x) > Band(x-3)*Thr */
            delta_energy =
                L_sub(Mpy_32_16(d2_fx[iBand - brickwall_dist], BW_thresh_brickwall[counter + 1]), d2_fx[iBand]);
            if (delta_energy > 0)
            {
                brickwall = 1; move16();
            }
            IF (brickwall)
            {
                BREAK;
            }
        }
        if (brickwall == 0)
        {
            *bw_idx = fs_idx;
        }
    }

    Dyn_Mem_Deluxe_Out();
}

