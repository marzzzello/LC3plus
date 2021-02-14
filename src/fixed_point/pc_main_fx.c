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


void processPCmain_fx(Word16 rframe, Word16 *bfi, Word16 prev_bfi, Word16 yLen, Word16 frame_dms, Word16 q_old_res_fx[],
                      Word16 *q_old_res_fx_exp, Word16 q_res_fx[], Word16 q_old_d_fx[], Word16 spec_inv_idx,
                      Word16 pitch_present, Word16 stab_fac, Word32 q_d_fx[], Word16 *q_fx_exp,
                      Word16 gg_idx, Word16 gg_idx_off, Word16 *prev_gg, Word16 *prev_gg_e, Word16 *BW_cutoff_idx_nf,
                      Word16 *prev_BW_cutoff_idx_nf, Word16 fac_ns_idx, Word16 *prev_fac_ns_fx, Word16 *pc_nbLostFramesInRow)
{
    Dyn_Mem_Deluxe_In(
        Word16 fac, fac_e;
    );
    
    fac = 32767; fac_e = 0;

    IF (sub(*bfi, 2) == 0)
    {
        processPCclassify_fx(pitch_present, frame_dms, q_old_d_fx, q_old_res_fx, yLen, spec_inv_idx, stab_fac, prev_bfi, bfi);
    }

    IF (sub(*bfi, 2) == 0)
    {
        processPCapply_fx(yLen, q_old_res_fx, q_old_res_fx_exp, q_res_fx, q_old_d_fx, spec_inv_idx,
                          &fac, &fac_e, q_d_fx, q_fx_exp, gg_idx, gg_idx_off, *prev_gg, *prev_gg_e, pc_nbLostFramesInRow);
    }

    IF (sub(*bfi, 1) != 0)
    {
        processPCupdate_fx(*bfi, yLen, q_old_res_fx, q_old_res_fx_exp, q_res_fx, spec_inv_idx, gg_idx, gg_idx_off, prev_gg,
                           prev_gg_e, rframe, BW_cutoff_idx_nf, prev_BW_cutoff_idx_nf, fac_ns_idx, prev_fac_ns_fx, fac, fac_e);
    }

    if (sub(*bfi, 2) != 0)
    {
        *pc_nbLostFramesInRow = 0;  move16();
    }

    Dyn_Mem_Deluxe_Out();
}


