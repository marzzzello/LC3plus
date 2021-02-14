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


void processPLCapply_fx(Word16 concealMethod, Word16 nbLostFramesInRow, Word16 bfi, Word16 prev_bfi,
                        Word16 frame_length, Word16 la_zeroes, const Word16 w[], Word16 x_fx[], Word16 ola_mem[],
                        Word16 *ola_mem_exp, Word16 q_old_d_fx[], Word16 *q_old_fx_exp, Word32 q_d_fx[],
                        Word16 *q_fx_exp, Word16 yLen, Word16 fs_idx, Word16 *damping, Word16 old_pitch_int,
                        Word16 old_pitch_fr, Word16 *ns_cum_alpha, Word16 *ns_seed, Word16 frame_dms, AplcSetup *plcAd,
                        Word8 *scratchBuffer)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
        Word32 *d2_fx;
        Word32 *q_old_d_fx32;
        Word32 *r_fx;
        Word32 *tdc_A_32;
        Word16  d2_fx_exp;
        Word16  r_fx_exp;
        Word16  Q_syn;
        Word8 * buffer_perBandEnergy, *buffer_preEmphasis, *buffer_InverseODFT, *buffer_Levinson, *buffer_tdc,
            *buffer_tdac, *buffer_phecu;
        Word16        y_e;             /*exponent of L_ecu_rec */
        Word16        tmp_is_trans[2]; /* may be  changed to a single variable */
        Word16        env_stab;
        
        Word16        n_bands, prev_bfi_plc2;
        const Word16 *band_offsets;
        Word32 *      L_ecu_rec; /*  local xtda  output is MAX_LEN -> input  buffer,
                                    as  tmp buffer for w32 fft MAX_LPROT */
    );

    d2_fx        = (Word32 *)scratchAlign(scratchBuffer, 0); /* Size = 4 * MAX_BANDS_NUMBER_PLC */
    q_old_d_fx32 = (Word32 *)scratchAlign(d2_fx, sizeof(*d2_fx) * MAX_BANDS_NUMBER_PLC); /* Size = 4 * MAX_BW */
    r_fx         = (Word32 *)scratchAlign(d2_fx, sizeof(*d2_fx) * MAX_BANDS_NUMBER_PLC); /* Size = 4 * (M + 1) */
    tdc_A_32     = (Word32 *)scratchAlign(r_fx, sizeof(*r_fx) * (M + 1));                /* Size = 4 * (M + 1) */

    L_ecu_rec = (Word32 *)scratchAlign(tdc_A_32, sizeof(*tdc_A_32) * (M + 1)); /* Size = 4 * MAX_LPROT bytes */

    buffer_perBandEnergy =
        (Word8 *)scratchAlign(q_old_d_fx32, sizeof(*q_old_d_fx32) * (MAX_BW)); /* Size = 2 * MAX_BANDS_NUMBER_PLC */
    buffer_preEmphasis =
        (Word8 *)scratchAlign(tdc_A_32, sizeof(*tdc_A_32) * (M + 1)); /* Size = 2 * MAX_BANDS_NUMBER_PLC */
    buffer_InverseODFT = buffer_preEmphasis;                          /* Size = 640 bytes */
    buffer_Levinson    = buffer_preEmphasis;                          /* Size = 4 * (M + 1) */
    buffer_tdc = scratchBuffer; /* Size = 2 * (MAX_PITCH + MAX_LEN/2 + MAX_LEN + MDCT_MEM_LEN_MAX + M + MAX_PITCH +
                                   MAX_LEN/2 + M + 1) bytes */

    buffer_tdac  = scratchBuffer; /* Size = 2 * MAX_LEN bytes */
    buffer_phecu = scratchBuffer; /* Size = 2 * MAX_LGW + 8 * MAX_LPROT + 12 * MAX_L_FRAME */
    /* Buffers overlap since they are not used at once */

    
#ifdef NONBE_PLC4_ADAP_DAMP
    UNUSED(ns_cum_alpha);
    UNUSED(ns_seed);
#endif

    /* Apply/Prepare PLC in bfi-case */
    IF (sub(bfi, 1) == 0)
    {
        SWITCH (concealMethod)
        {
        case 2:
            ASSERT(frame_dms == 100);
            /* call phaseEcu */
            env_stab        = 32767; move16();                 /* 1.0=stable , 0.0=dynamic Q15*/
            tmp_is_trans[0] = plcAd->PhECU_short_flag_prev; move16();
            tmp_is_trans[1] = plcAd->PhECU_short_flag_prev; move16();
            
            ASSERT(prev_bfi == 0 || prev_bfi == 1|| prev_bfi == 2);  /*PC prev_bfi has three states */
            prev_bfi_plc2 = prev_bfi; move16();
            if (sub(prev_bfi_plc2, 2) == 0) 
            {
               prev_bfi_plc2 = 0; move16();
            }  

            ASSERT(prev_bfi_plc2 == 0 || prev_bfi_plc2 == 1); /*PhEcu does not accept prev_bfi == 2 */
            IF(prev_bfi_plc2 == 0)
            { /* convert pitch lag info at current fs to a normalized fractional bin-frequency   */
               plcAd->PhECU_f0hzLtpBinQ7 = plc_phEcuSetF0Hz_fx(fs_idx, old_pitch_int, old_pitch_fr);  move16();







               /* first bfi frame calc decoded pcm  energy 16,16 ms, in 26 ms buffer separated by 10 ms*/
                
               {   /* compute energy normalization needed for concealment method 2  Xavg  and transient analysis */
                   
                    
                   /* left   */                
                  processPLCUpdateXFP_w_E_hist_fx(0, 0,
                     &(plcAd->x_old_tot_fx[ sub(plcAd->max_len_pcm_plc , add(num_FsByResQ0[fs_idx],rectLengthTab[fs_idx]  )) ]), plcAd->q_fx_old_exp,0,  
                  
                     fs_idx,
                     &plcAd->PhECU_L_oold_xfp_w_E_fx, &plcAd->PhECU_oold_xfp_w_E_exp_fx,
                     &plcAd->PhECU_L_old_xfp_w_E_fx, &plcAd->PhECU_old_xfp_w_E_exp_fx,
                     &plcAd->PhECU_oold_Ltot_exp_fx, &plcAd->PhECU_old_Ltot_exp_fx);

                  /* right  */
                  processPLCUpdateXFP_w_E_hist_fx(0, 0, plcAd->PhECU_xfp_fx, plcAd->PhECU_xfp_exp_fx, 
                     plcAd->PhECU_margin_xfp, fs_idx,
                     &plcAd->PhECU_L_oold_xfp_w_E_fx, &plcAd->PhECU_oold_xfp_w_E_exp_fx,
                     &plcAd->PhECU_L_old_xfp_w_E_fx, &plcAd->PhECU_old_xfp_w_E_exp_fx,
                     &plcAd->PhECU_oold_Ltot_exp_fx, &plcAd->PhECU_old_Ltot_exp_fx);               

               }
            }

            hq_phase_ecu_fx(
                plcAd->PhECU_xfp_fx,       /* i :  only valid first Bfi frame , buffer of previous synt signal length */
                L_ecu_rec,                 /* o  : reconstructed frame in folded tda domain xtda  Word32  Q x     */
                &plcAd->PhECU_time_offs,   /* i/o: Sample offset for consecutive frame losses*/
                plcAd->PhECU_X_sav_fx,     /* i(prev_bfi==1)/o(prev_bfi==0): Stored Complex spectrum of prototype frame */
                &plcAd->PhECU_X_savQ_fx,   /* i/o: Q value of stored spectrum                */
                &plcAd->PhECU_num_plocs,   /* i/o: Number of identified peaks                */
                plcAd->PhECU_plocs,        /* i/o: Peak locations   Q0                         */
                plcAd->PhECU_f0est,        /* i/o: Interpolated peak locations           Q16 */
                env_stab,                  /* i  : Envelope stability parameter              */
                plcAd->PhECU_f0hzLtpBinQ7, /* i:  LTP bin frequency in normalized Hz  Q7 */
                plcAd->norm_corrQ15_fx,    /* i : correlation for lag at f0hzLtpBinQ7 */
                prev_bfi_plc2,                  /* i   : indicating burst frame error             */
                tmp_is_trans,              /* i   : flags indicating previous transient frames */
                plcAd->PhECU_mag_chg_1st,  /* i/o: per band magnitude modifier for transients*/
                NULL,                      /*  o: dbg per band magnitude modifier, incl. burst attenuation   */
                plcAd->PhECU_Xavg,         /* i/o: Frequency group average gain to fade to   */
                &plcAd->PhECU_beta_mute,   /* o   : Factor for long-term mute                */
                fs_idx,                    /* i  : Encoded bandwidth   "nb(0),WB,sWB,WB,FB"  */
                frame_length,              /* i   : frame length                             */
                NULL ,                      /* o  :  seed synch dbg                          */ 
                NULL ,                      /* o  :  evolved Spectrum   dbg                  */ 
                plcAd->PhECU_t_adv,       /* i  : time adjustment excluding time_offs         */
                 PhECU_wins[fs_idx][2], /* i:  2 ms initial part pre_tda = mdct-ana */
                 PhECU_wins[fs_idx][1], /* i:  16 ms pretda combined part  IWHR+MDCT-ana  */
                PhECU_wins[fs_idx][0], 
                plcAd->PhECU_xfp_exp_fx, 
                plcAd->max_lprot, 
                plcAd->max_plocs,
                plcAd->PhECU_L_oold_xfp_w_E_fx,plcAd->PhECU_oold_xfp_w_E_exp_fx, plcAd->PhECU_oold_Ltot_exp_fx,
                plcAd->PhECU_oold_grp_shape_fx,
                plcAd->PhECU_L_old_xfp_w_E_fx,plcAd->PhECU_old_xfp_w_E_exp_fx, plcAd->PhECU_old_Ltot_exp_fx,
                plcAd->PhECU_old_grp_shape_fx,
                plcAd->PhECU_margin_xfp,
                buffer_phecu);
 
            y_e = 18;  move16();  /*  the  fixed exponent (exp)  from Lecu_rec  from PhaseECU is 18    */

            Processing_ITDA_WIN_OLA(
                L_ecu_rec,   /* i:     X_TDA buffer data   =  "y"  DCT-IV output */
                &y_e,        /* i/o:    x_tda exponent  "y_e"  */
                w,           /* i:     window coefficients including normalization of sqrt(2/N) and scaled by 2^4 */
                ola_mem,     /* i/o:  overlap add memory */
                ola_mem_exp, /* i/o:  overlap add exponent */
                x_fx,        /* o:    time signal out */
                LowDelayShapes_n960_len[fs_idx],                         /* i:   window length */
                frame_length,                                            /* i:   block size */
                sub(frame_length, LowDelayShapes_n960_la_zeroes[fs_idx]) /* i:   overlap add buffer size */
            ); 
           *q_fx_exp = y_e;  move16();   /*  assign updated Q */


            BREAK;

        case 3:
            IF (sub(nbLostFramesInRow, 1) == 0)
            {
                plcAd->tdc_fract = old_pitch_fr; move16();
                n_bands          = s_min(frame_length, MAX_BANDS_NUMBER_PLC);
                SWITCH (frame_dms)
                {
                case 25:
                    band_offsets = bands_offset_lin_2_5ms[fs_idx];  move16();
                    IF (sub(fs_idx, 4) == 0)
                    {
                        n_bands = 60;  move16();
                    }
                    BREAK;
                case 50:
                    band_offsets = bands_offset_lin_5ms[fs_idx]; move16();
                    IF (sub(fs_idx, 2) == 0)
                    {
                        n_bands = 40; move16();
                    }
                    BREAK;
                case 100:
                    band_offsets = bands_offset_lin[fs_idx]; move16();
                    BREAK;
                }

                FOR (i = 0; i < yLen; i++)
                {
                    q_old_d_fx32[i] = L_deposit_h(q_old_d_fx[i]);
                }

                /* LPC Analysis */
                /* calculate per band energy*/
                processPerBandEnergy_fx(d2_fx, &d2_fx_exp, q_old_d_fx32, *q_old_fx_exp, band_offsets, fs_idx, n_bands,
                                        1, frame_dms, buffer_perBandEnergy);

                /* calculate pre-emphasis */
                processPreEmphasis_fx(d2_fx, &d2_fx_exp, fs_idx, n_bands, frame_dms, buffer_preEmphasis);

                /* inverse ODFT */
                processInverseODFT_fx(r_fx, &r_fx_exp, d2_fx, d2_fx_exp, n_bands, plcAd->tdc_lpc_order, buffer_InverseODFT);

                /* lag windowing */
                processLagwin_fx(r_fx, lag_win[fs_idx], plcAd->tdc_lpc_order);

                /* Levinson Durbin */
                processLevinson_fx(tdc_A_32, r_fx, plcAd->tdc_lpc_order, NULL, NULL, buffer_Levinson);

                /* 32Q27 -> 16Qx */
                processPLCLpcScaling_fx(tdc_A_32, plcAd->tdc_A, add(plcAd->tdc_lpc_order, 1));
            }

            /* call TD-PLC */
            processTimeDomainConcealment_Apply_fx(
                old_pitch_int, plcAd->tdc_preemph_fac, plcAd->tdc_A, plcAd->tdc_lpc_order, plcAd->x_old_tot_fx, frame_length, frame_dms,
                fs_idx, nbLostFramesInRow, sub(frame_length, la_zeroes), plcAd->stab_fac, &plcAd->tdc_fract,
                &plcAd->tdc_seed, &plcAd->tdc_gain_p, &plcAd->tdc_gain_c, &plcAd->tdc_cum_damp, x_fx, &Q_syn, damping,
                plcAd->max_len_pcm_plc, buffer_tdc);

            /* exponent of TD-PLC output */
            Q_syn     = add(Q_syn, sub(15, plcAd->q_fx_old_exp));
            *q_fx_exp = sub(15, Q_syn); move16();

            /* TDAC */
            processTdac_fx(ola_mem, ola_mem_exp, x_fx, *q_fx_exp, w, la_zeroes, frame_length, buffer_tdac);
            BREAK;

        case 4:
            *q_fx_exp = *q_old_fx_exp; move16();

            /* call Noise Substitution */
#ifndef NONBE_PLC4_ADAP_DAMP
            processPLCNoiseSubstitution_fx(q_d_fx, q_old_d_fx, yLen, nbLostFramesInRow, plcAd->stab_fac, frame_dms,
                                           damping, ns_cum_alpha, ns_seed);
#else
            processPLCNoiseSubstitution_fx(q_d_fx, q_old_d_fx, yLen);
#endif
            BREAK;

        default: ASSERT(!"Unsupported PLC method!");
        } /* switch (converalMethod)*/
    }     /* if (bfi) */

    Dyn_Mem_Deluxe_Out();
}


