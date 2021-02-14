/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"

#include "functions.h"


static Word16 spectral_centroid_fx_lc(Word16 old_scf_q[], const Word16 *band_offsets, Word16 frame_length,
                                      Word16 fs_idx, Word8 *scratchBuffer);

void processPLCclassify_fx(Word16 plcMeth, Word16 *concealMethod, Word16 *nbLostFramesInRow, Word16 bfi,
                           Word16 ltpf_mem_pitch_int, Word16 frame_length, Word16 frame_dms, Word16 fs_idx, Word16 yLen,
                           Word16 q_old_d_fx[], const Word16 *band_offsets, AplcSetup *plcAd, Word8 *scratchBuffer)
{
    Dyn_Mem_Deluxe_In(
        Word16 scQ15;
        Word32 class;
    );

    BASOP_sub_sub_start("PLC::processPLCclassify_fx");

    UNUSED(yLen);
    UNUSED(q_old_d_fx);

    if (plcAd)
    {
        plcAd->norm_corrQ15_fx = 0; move16();
    }
    IF (sub(bfi, 1) == 0)
    {
        /* increase counter of lost-frames-in-a-row */
        *nbLostFramesInRow = add(*nbLostFramesInRow, 1);

        IF (sub(*nbLostFramesInRow, 1) == 0)
        {
            *concealMethod = plcMeth; move16();

            IF (sub(plcMeth, 1) == 0)
            {
                IF (ltpf_mem_pitch_int > 0)
                {
                    *concealMethod = 3; move16(); /* TD-PLC */
                    /* no classifier features needed for 5ms mode, as PhaseECU-5ms is disabled */
                    IF (sub(frame_dms, 100) == 0)
                    {
                        /* Calculate Features */
                        plcAd->norm_corrQ15_fx = plc_xcorr_lc_fx(plcAd->x_old_tot_fx, plcAd->max_len_pcm_plc,
                                                                 ltpf_mem_pitch_int, frame_length, fs_idx);
                        scQ15 = spectral_centroid_fx_lc(plcAd->old_scf_q, band_offsets, frame_length, fs_idx,
                                                        scratchBuffer);

                        /* Classify */
                        class = L_mult(plcAd->norm_corrQ15_fx, 7640);
                        class = L_mac(class, scQ15, -32768);
                        class = L_add_sat(class, -335020208);

                        IF (class <= 0)
                        {
                            *concealMethod = 2; move16(); /* Phase ECU selected */
                        }
                    }
                }
                ELSE
                {
                    *concealMethod = 4; move16(); /* Noise Substitution */
                }
            }
        }
    }

    Dyn_Mem_Deluxe_Out();
    BASOP_sub_sub_end();
}


Word16 spectral_centroid_fx_lc(Word16 old_scf_q[], const Word16 *band_offsets, Word16 frame_length, Word16 fs_idx,
                               Word8 *scratchBuffer)
{
    Dyn_Mem_Deluxe_In(
        Counter i, j;
        Word32  den32, num32, tmp32;
        Word16  s, sc, fac, freq, inv, startfreq, stopfreq;
        Word16 *old_scf_q_mod;
        Word16 *old_scf_q_mod_exp;
    );
    BASOP_sub_sub_start("PLC::spectral_centroid_fx_lc");

    old_scf_q_mod     = (Word16 *)scratchAlign(scratchBuffer, 0);                          /* Size = 2 * M */
    old_scf_q_mod_exp = (Word16 *)scratchAlign(old_scf_q_mod, sizeof(*old_scf_q_mod) * M); /* Size = 2 * M */

    /* Linear Domain */
    FOR (i = 0; i < M; i++)
    {
        old_scf_q_mod[i] = BASOP_Util_InvLog2_16(old_scf_q[i], &old_scf_q_mod_exp[i]);
    }

    /* De-emphasis */
    FOR (i = 0; i < M; i++)
    {
        old_scf_q_mod[i]     = mult(old_scf_q_mod[i], lpc_warp_dee_emphasis[fs_idx][i]);      move16();
        old_scf_q_mod_exp[i] = add(old_scf_q_mod_exp[i], lpc_warp_dee_emphasis_e[fs_idx][i]); move16();
    }

    den32 = 1; move16();
    num32 = 0; move16();
    inv   = div_s(1, frame_length);

    FOR (i = 0; i < M; i++)
    {
        freq      = 0; move16();
        startfreq = add(band_offsets[i * 4], 1);
        stopfreq  = band_offsets[i * 4 + 4];
        FOR (j = startfreq; j <= stopfreq; j++)
        {
            freq = add(freq, j);
        }

        tmp32 = L_mult(inv, freq);
        s     = norm_l(tmp32);
        tmp32 = L_mult(old_scf_q_mod[i], extract_h(L_shl(tmp32, s)));

        num32 = L_add(num32, L_shl(tmp32, add(add(-15, old_scf_q_mod_exp[i]), sub(15, s))));
        den32 = L_add(den32, L_shl(L_mult(old_scf_q_mod[i], stopfreq - startfreq + 1), old_scf_q_mod_exp[i]));
    }

    s = norm_l(den32);
    s = sub(16, s);

    sc = div_s(extract_l(L_shr(num32, s)), extract_l(L_shr(den32, s)));

    SWITCH (fs_idx)
    {
    case 0:
        fac = 5461; move16();
        BREAK;
    case 1:
        fac = 10922; move16();
        BREAK;
    case 2:
        fac = 16384; move16();
        BREAK;
    case 3:
        fac = 21845; move16();
        BREAK;
    default:         /* case 4: */
        fac = 32767; move16();
        BREAK;
    }
    sc = round_fx(L_mult(sc, fac));

    Dyn_Mem_Deluxe_Out();
    BASOP_sub_sub_end();
    return sc;
}


