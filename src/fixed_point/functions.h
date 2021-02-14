/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "basop_util.h"
#include "basop_mpy.h"
#include "constants.h"
#include "defines.h"
#include "dynmem.h"
#include "lc3.h"
#include "rom_basop_util.h"
#include "setup_dec_lc3.h" /* for decoder state handle ptr  */
#include "setup_enc_lc3.h" /* for encoder state handle ptr  */
#include "stl.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PRINTF printf /* C console debug,     */
#define ASSERT(test) assert(test)

void WarnMsg(char *msg);
void ExitMsg(char *msg);
void AssertMsg(int true_flag, char *msg);

void processTnsDecoder_fx(Word16 rc_idx[], Word32 x[], Word16 xLen, Word16 order[], Word16 *x_e, Word16 BW_stopband_idx,
                          Word16 frame_dms, Word8 *scratchBuffer);

void processTnsCoder_fx(Word16 *bits, Word16 indexes[], Word32 x[], Word16 BW_cutoff_idx, Word16 order[],
                        Word16 *numfilters, Word16 enable_lpc_weighting, Word16 nSubdivisions, Word16 frame_dms,
                        Word16 max_len, Word8 *scratchBuffer);

void TnsUpdate_fx(Word16 *tns_numfilters, Word16 *long_startfreq, const Word16 **subdiv_startfreq,
                  const Word16 **subdiv_stopfreq, Word16 frame_length, Word16 tab_idx);

void processResidualDecoding_fx(Word32 spectrum[], Word16 spectrum_e, Word16 L_spec, Word16 prm[], Word16 resQBits);

void processPreemph_fx(Word16 *x, Word16 len, Word16 mem[], Word16 memLen, Word16 out[], Word16 *outLen, Word16 mu);

void processNoiseFilling_fx(Word32 xq[], Word16 nfseed, Word16 xq_e, Word16 fac_ns_idx, Word16 BW_cutoff_idx,
                            Word16 frame_dms, Word16 fac_ns_pc, Word16 spec_inv_idx, Word8 *scratchBuffer);

#ifdef NONBE_LOW_BR_NF_TUNING
void processNoiseFactor_fx(Word16 *fac_ns_idx, Word16 x_e, Word32 x[], Word16 xq[], Word16 gg, Word16 gg_e,
                           Word16 BW_cutoff_idx, Word16 frame_dms, Word16 target_bytes, Word8 *scratchBuffer);
#else
void processNoiseFactor_fx(Word16 *fac_ns_idx, Word16 x_e, Word32 x[], Word16 xq[], Word16 gg, Word16 gg_e,
                           Word16 BW_cutoff_idx, Word16 frame_dms, Word8 *scratchBuffer);
#endif

void processMdctShaping_fx(Word32 x[], Word16 scf[], Word16 scf_exp[], const Word16 bands_offset[], Word16 fdns_npts);

void processScfScaling(Word16 scf_exp[], Word16 fdns_npts, Word16 *x_e);

void processMdct_fx(Word16 x[], Word16 x_exp, Word16 N, const Word16 w[], Word16 wLen, Word16 mem[], Word16 memLen,
                    Word32 y[], Word16 *y_e, Word8 *scratchBuffer);

void processLevinson_fx(Word32 *lpc, Word32 *ac, Word16 N, Word16 *rc, Word32 *pred_err, Word8 *scratchBuffer);

void lpc2rc(Word32 *lpc, Word16 *rc, Word16 N);

void ProcessingIMDCT(Word32 y[], Word16 *y_e, const Word16 w[], Word16 mem[], Word16 *mem_e, Word16 x[], Word16 wLen,
                     Word16 N, Word16 memLen, Word16 frame_dms,
                     Word16 concealMethod, Word16 bfi, Word16 prev_bfi, Word16 nbLostFramesInRow, AplcSetup *plcAd,
                     Word8 *scratchBuffer);

void dct_IV(Word32 *pDat, Word16 *pDat_e, Word16 L, Word32 *workBuffer);

Word16 find_last_nz_pair(const Word16 x[], Word16 length);

PvqEntry_fx mpvq_index_fx(const Word16 *vec_in, Word16 dim_in, Word16 k_val_local);
void        mpvq_deindex_fx(const PvqEntry_fx *entry, UWord32 *h_mem, Word16 *vec_out);

void pvq_enc_search_fx(const Word16 *x,       /* i:   target vector to quantize             Qin      */
                       Word16 *      rt_far,  /* o:   outl_far o, raw pulses  (non-scaled short) Q0      */
                       Word16 *      rt_near, /* o:   outl_near o, raw pulses  (non-scaled short) Q0      */
                       Word16 *      rtA,     /* o:   A section  raw pulses  (non-scaled short) Q0    */
                       Word16 *      rtB,     /* o:   B section  raw pulses  (non-scaled short) Q0     */
                       Word32 *L_corr,      /* o:   4 un-normalized correlation sums for outl_far, near, outl, A, AB  */
                       Word32 *L_search_en, /* o:   4 energy sums for out_far, outl_near, A, AB  */
                       Word16 *pulses_fin,  /* i:   number of allocated pulses  to outl A, AB section  */
                       Word16 *pulses_proj, /* i:   number of proj pulses  pulses to outl, A , AB section     */

                       const Word16 dim, /* i:   Length of vector */
                       const Word16 dimA /* i:   Length of vector A section */
);

PvqEntry_fx get_size_mpvq_calc_offset_fx(Word16 dim_in, Word16 k_val_in, UWord32 *h_mem);

Word16 pvq_dec_deidx_fx(                          /* out BER detected 1 , ok==0 */
                        Word16 *      y,          /* o:   decoded vector (non-scaled int)  */
                        const Word16  k_val,      /* i:   number of allocated pulses       */
                        const Word16  dim,        /* i:   Length of vector                 */
                        const Word16  LS_ind,     /* i; lS index              1 bit      */
                        const UWord32 UL_MPVQ_ind /* i; MPVQ  index                      */
);

void pvq_dec_en1_norm_fx(                            /*  */
                         Word16 *      xq,           /* o:   normalized decoded vector (Q15)             */
                         const Word16 *y,            /* i:   decoded vector (non-scaled int)  */
                         const Word16  kval_max,     /* i:   max possible K   in Q0 kO or kA   */
                         const Word16  dim,          /* i:   Length of vector                 */
                         const Word16  neg_glob_gain /* i:   a Global Gain   (negated to fit 1.0 in Q15 as -1.0) */
);

void pvq_dec_en1_normQ14_fx(                         /*  Have to be used by both encoder and decoder */
                            Word16 *      xq,        /* o:   en1 normalized decoded vector (Q14)        */
                            const Word16 *y,         /* i:   decoded vector (non-scaled int)  */
                            const Word16  k_val_max, /* i:   max possible K   in Q0 kO or kA   */
                            const Word16  dim        /* i:   Length of vector                 */
);

void pvq_dec_scale_vec_fx(const Word16 *inQ14, Word16 adjGainQ13, Word16 *vecQ14);

Word16 processAriEncoder_fx(UWord8 *bytes, Word16 bp_side, Word16 mask_side, Word16 nbbits, Word16 xq[],
                            Word16 *tns_order, Word16 tns_numfilters, Word16 *tns_idx, Word16 lastnz,
                            Word16 *codingdata, Word8 *resBits, Word16 numResBits, Word16 lsbMode,
                            Word16 enable_lpc_weighting, Word8 *scratchBuffer);

void processAriDecoder_fx(UWord8 *bytes, Word16 *bp_side, Word16 *mask_side, Word16 nbbits, Word16 L_spec,
                          Word16 fs_idx, Word16 enable_lpc_weighting, Word16 tns_numfilters, Word16 lsbMode,
                          Word16 lastnz, Word16 *bfi, Word16 *tns_order, Word16 fac_ns_idx, Word16 gg_idx,
                          Word16 frame_dms,
                          Word16 n_pc, Word16 be_bp_left, Word16 be_bp_right, Word16 enc, Word16 *spec_inv_idx, Word16 *b_left,
                          Word16 *resBits, Word16 *x, Word16 *nf_seed, Word16 *resQdata, Word16 *tns_idx,
                          Word16 *zero_frame, Word8 *scratchBuffer);

void processAriDecoderScaling_fx(Word16 *data16, Word16 dataLen, Word32 *data32, Word16 *data_e);

void processApplyGlobalGain_fx(Word32 x[], Word16 *x_e, Word16 xLen, Word16 global_gain_idx, Word16 global_gain_off);

void processPerBandEnergy_fx(Word32 *d2_fx, Word16 *d2_fx_exp, Word32 *d_fx, Word16 d_fx_exp,
                             const Word16 *band_offsets, Word16 fs_idx, Word16 n_bands, Word16 linear, Word16 frame_dms,
                             Word8 *scratchBuffer);

void processDetectCutoffWarped_fx(Word16 *bw_idx, Word32 *d2_fx, Word16 d2_fx_exp, Word16 fs_idx, Word16 frame_dms);

void process_resamp12k8_fx(Word16 x[], Word16 x_len, Word16 mem_in[], Word16 mem_in_len, Word32 mem_50[],
                           Word16 mem_out[], Word16 mem_out_len, Word16 y[], Word16 *y_len, Word16 fs_idx,
                           Word16 frame_dms, Word8 *scratchBuffer);

void process_olpa_fx(Word16 *mem_s6k4_exp, Word16 mem_s12k8[], Word16 mem_s6k4[], Word16 *pitch, Word16 *s12k8,
                     Word16 len, Word16 *normcorr, Word16 *mem_pitch, Word16 s12k8_exp, Word8 *scratchBuffer);

void process_ltpf_coder_fx(Word16 *bits, Word16 ol_pitch, Word16 ltpf_enable, Word16 *old_wsp_exp, Word16 *old_wsp,
                           Word16 old_wsplen, Word16 *param, Word16 *wsp, Word16 len, Word16 *mem_normcorr,
                           Word16 *mem_mem_normcorr, Word16 ol_normcorr, Word16 *mem_ltpf_on, Word16 *mem_ltpf_pitch,
                           Word16 wsp_exp, Word16 frame_dms, Word8 *scratchBuffer);

void process_ltpf_decoder_fx(Word16 *x_e, Word16 L_frame, Word16 old_x_len, Word16 fs_idx, Word16 old_y_len,
                             Word16 *old_e, Word16 *x, Word16 *old_x, Word16 *y, Word16 *old_y, Word16 ltpf,
                             Word16 ltpf_active, Word16 pitch_index, Word16 *old_pitch_int, Word16 *old_pitch_fr,
                             Word16 *old_gain, Word16 *mem_ltpf_active, Word16 scale_fac_idx, Word16 bfi,
                             Word16 concealMethod,
                             Word16 damping, Word16 *old_scale_fac_idx, Word8 *scratchBuffer);

void attack_detector_fx(LC3_Enc *enc, EncSetup *setup, Word16 *input, Word16 input_scaling, void *scratch);

void processSnsComputeScf_fx(Word32 *d2_fx, Word16 d2_fx_exp, Word16 fs_idx, Word16 n_bands, Word16 *scf,
                             Word16 scf_smoothing_enabled, Word8 *scratchBuffer);

void processSnsQuantizeScfEncoder_fx(Word16  scf[],     /* i: input scf M */
                                     Word32 *L_prm_idx, /* o: indeces . negative == unused */
                                     Word16 *scf_q,     /* o: quantized scf M */
                                     Word8 * scratchBuffer);

Word16 processSnsQuantizeScfDecoder_fx(                                       /* o: BER flag */
                                       Word32 *L_prm_idx,                     /* i: indeces */
                                       Word16 scf_q[], Word8 *scratchBuffer); /* o:  M */

void processSnsInterpolateScf_fx(Word16 *scf_q, Word16 mdct_scf[], Word16 mdct_scf_exp[], Word16 inv_scf,
                                 Word16 n_bands, Word8 *scratchBuffer);

void processPLCmain_fx(Word16 plcMeth, Word16 *concealMethod, Word16 *nbLostFramesInRow, Word16 bfi, Word16 prev_bfi,
                       Word16 frame_length, Word16 la_zeroes, const Word16 w[], Word16 x_fx[], Word16 ola_mem[],
                       Word16 *ola_mem_exp, Word16 q_old_d_fx[], Word16 *q_old_fx_exp, Word32 q_d_fx[],
                       Word16 *q_fx_exp, Word16 yLen, Word16 fs_idx, const Word16 *band_offsets, Word16 *damping,
                       Word16 old_pitch_int, Word16 old_pitch_fr, Word16 *ns_cum_alpha, Word16 *ns_seed,
                       AplcSetup *plcAd, Word16 frame_dms, Word8 *scratchBuffer);

void processPLCupdate_fx(AplcSetup *plcAd, Word16 x_fx[], Word16 q_fx_exp, Word16 concealMethod, Word16 frame_length,
                         Word16 fs_idx, Word16 *nbLostFramesInRow, Word16 *prev_prev_bfi, Word16 *prev_bfi, Word16 bfi,
                         Word16 scf_q[], Word16 ola_mem_fx[], Word16 ola_mem_fx_exp, Word16 *ns_cum_alpha);

void processPLCupdateSpec_fx(Word16 q_old_d_fx[], Word16 *q_old_fx_exp, Word32 q_d_fx[], Word16 *q_fx_exp, Word16 yLen);


void processPLCspec2shape_fx(Word16 prev_bfi, Word16 bfi, Word16 q_old_d_fx[], Word16 yLen, Word16 *PhECU_oold_grp_shape_fx, Word16 *PhECU_old_grp_shape_fx);


Word32 winEnCalc(const Word16 *, const Word16 , const Word16 *, const Word16, const Word16,Word16 *);

void processPLCUpdateAfterIMDCT_fx(Word16 x_fx[], Word16 q_fx_exp, Word16 concealMethod, Word16 xLen, Word16 fs_idx,
                                   Word16 *nbLostFramesInRow, Word16 *prev_prev_bfi, Word16 *prev_bfi, Word16 bfi, Word16 scf_q[],
                                   Word16 *ns_cum_alpha, AplcSetup *plcAd);

void processPLCclassify_fx(Word16 plcMeth, Word16 *concealMethod, Word16 *nbLostFramesInRow, Word16 bfi,
                           Word16 ltpf_mem_pitch_int, Word16 frame_length, Word16 frame_dms, Word16 fs_idx, Word16 yLen,
                           Word16 q_old_d_fx[], const Word16 *band_offsets, AplcSetup *plcAd, Word8 *scratchBuffer);

void processPLCapply_fx(Word16 concealMethod, Word16 nbLostFramesInRow, Word16 bfi, Word16 prev_bfi,
                        Word16 frame_length, Word16 la_zeroes, const Word16 w[], Word16 x_fx[], Word16 ola_mem[],
                        Word16 *ola_mem_exp, Word16 q_old_d_fx[], Word16 *q_old_fx_exp, Word32 q_d_fx[],
                        Word16 *q_fx_exp, Word16 yLen, Word16 fs_idx, Word16 *damping, Word16 old_pitch_int,
                        Word16 old_pitch_fr, Word16 *ns_cum_alpha, Word16 *ns_seed, Word16 frame_dms, AplcSetup *plcAd,
                        Word8 *scratchBuffer);

#ifndef NONBE_PLC4_ADAP_DAMP
void processPLCNoiseSubstitution_fx(Word32 spec[], Word16 spec_prev[], Word16 L_spec, Word16 nbLostFramesInRow,
                                    Word16 stabFac, Word16 frame_dms, Word16 *alpha, Word16 *cum_alpha, Word16 *seed);
#else
void processPLCNoiseSubstitution_fx(Word32 spec[], Word16 spec_prev[], Word16 L_spec);
void processPLCDampingScrambling_main_fx(Word16 bfi, Word16 concealMethod, Word16 ns_nbLostFramesInRow,
                                         Word16 pc_nbLostFramesInRow, Word16 *ns_seed, Word16 *pc_seed, Word16 pitch_present_bfi1,
                                         Word16 pitch_present_bfi2, Word32 spec[], Word16 *q_fx_exp, Word16 *q_old_d_fx,
                                         Word16 *q_old_fx_exp, Word16 L_spec, Word16 stabFac, Word16 frame_dms,
                                         Word16 *cum_fading_slow, Word16 *cum_fading_fast, Word16 *alpha, Word16 spec_inv_idx);
void processPLCDampingScrambling_fx(Word32 spec[], Word16 L_spec, Word16 nbLostFramesInRow, Word16 stabFac,
                                    Word16 pitch_present, Word16 frame_dms, Word16 *cum_fading_slow,
                                    Word16 *cum_fading_fast, Word16 *alpha, Word16 *seed, Word16 spec_inv_idx);
#endif

void processLagwin_fx(Word32 r[], const Word32 w[], Word16 m);

void processInverseODFT_fx(Word32 *r_fx, Word16 *r_fx_exp, Word32 *d2_fx, Word16 d2_fx_exp, Word16 n_bands,
                           Word16 lpc_order, Word8 *scratchBuffer);

void processPreEmphasis_fx(Word32 *d2_fx, Word16 *d2_fx_exp, Word16 fs_idx, Word16 n_bands, Word16 frame_dms, Word8 *scratchBuffer);

void processPLCLpcScaling_fx(Word32 tdc_A_32[], Word16 tdc_A_16[], Word16 m);

#ifdef BE_MOVED_STAB_FAC
void processPLCcomputeStabFac_main(Word16 scf_q[], Word16 old_scf_q[], Word16 old_old_scf_q[], Word16 bfi, Word16 prev_bfi,
                              Word16 prev_prev_bfi, Word16 *stab_fac);
void processPLCcomputeStabFac(Word16 scf_q[], Word16 old_scf_q[], Word16 prev_bfi, Word16 *stab_fac);
#endif

void processPLCUpdateXFP_w_E_hist_fx(Word16 prev_bfi, Word16 bfi, 
                                     Word16 *xfp_fx, Word16 xfp_exp_fx,Word16 margin_xfp, 
                                     Word16 fs_idx, 
                                     Word32 *L_oold_xfp_w_E_fx, Word16 *oold_xfp_w_E_exp_fx, 
                                     Word32 *L_old_xfp_w_E_fx, Word16 *old_xfp_w_E_exp_fx,                                   
                                     Word16 *oold_Ltot_exp_fx, Word16 *old_Ltot_exp_fx);    
void processTimeDomainConcealment_Apply_fx(const Word16 pitch_int, const Word16 preemphFac_fx, const Word16 *A_fx,
                                           const Word16 lpc_order, const Word16 *pcmbufHist_fx, const Word16 frame_length,
                                           const Word16 frame_dms, const Word16 fs_idx, const Word16 nbLostFramesInRow,
                                           const Word16 overlap, const Word16 stabFac_fx, Word16 *fract,
                                           Word16 *seed_fx, Word16 *gain_p_fx, Word32 *gain_c_fx, Word16 *cum_alpha,
                                           Word16 *synth_fx, Word16 *Q_syn, Word16 *alpha, Word16 max_len_pcm_plc,
                                           Word8 *scratchBuffer);

void processTdac_fx(Word16 *ola_mem, Word16 *ola_mem_exp, const Word16 *synth, const Word16 synth_exp,
                    const Word16 *win, const Word16 la_zeroes, const Word16 frame_len, Word8 *scratchBuffer);

void plc_phEcu_F0_refine_first_fx(Word16 *plocs, const Word16 n_plocs_in, Word32 *L_f0est,
                                  const Word16 stPhECU_f0hzLtpBinQ7, const Word16 stPhECU_f0gainLtpQ15,
                                  const Word16 nSubm);
void plc_phEcu_LF_peak_analysis_fx(Word16 *plocs, Word16 *n_plocs, Word32 *L_f0estQ16, const Word16 *mag,
                                   const Word16 stPhECU_f0hzLtpBinQ7, const Word16 stPhECU_f0gainLtpQ15,
                                   const Word16 nSubm, Word16 maxPlocs, Word8 *scratchBuffer);

Word16 plc_phEcuSetF0Hz_fx(Word16 fs_idx, Word16 old_pitch_int, Word16 old_pitch_fr);

void create_sin2_taper_fx(Word16 *, Word16, Word16);

void plc_phEcu_initWord16(Word16 *     vec,   /*i/o : vector pointer             */
                          const Word16 value, /*i   : short initialization value */
                          const Word16 len);  /*i   : number of elements         */

Word16 plc_phEcu_ratio_fx(const Word32, const Word32, Word16 *);

void plc_phEcu_minval_fx(const Word16 *inp,  /* i  : vector       */
                         const Word16  len,  /* i  : length       */
                         Word16 *      minvalPtr); /* o  : min  value Ptr    */

void plc_phEcu_maxval_fx(const Word16 *inp,  /* i  : vector     */
                         const Word16  len,  /* i  : length     */
                         Word16 *      maxvalPtr); /* o  : *maxvalPtr */

void Scale_sig_sat(Word16       x[],   /* i/o: signal to scale,  possibly saturated      Qx        */
                   const Word16 lg,    /* i  : size of x[]                     Q0        */
                   const Word16 exp0); /* i  : exponent: x = round(x << exp)   Qx ?exp  */

void Processing_ITDA_WIN_OLA(Word32 L_x_tda[], Word16 *y_e, const Word16 w[], Word16 mem[], Word16 *mem_e, Word16 x[],
                             Word16 wLen, Word16 N, Word16 memLen);

void trans_burst_ana_fx(const Word16 *xfp,     /* i  : Input signal                                       Qspec */
                        Word16 *      mag_chg, /* o  : Magnitude modification                             Q15 */
                        Word16 *ph_dith, /* o  : Phase dither, 2*PI is not included (Q15, i.e., between 0.0 and 1.0) */
                        Word16 *mag_chg_1st,          /* i/o: per band magnitude modifier for transients         Q15 */
                        const Word16 output_frame,    /* i  : Frame length                                           */
                        const Word16 time_offs,       /* i  : Time offset (integral multiple of output_frame)        */
                        const Word16 est_stab_content, /* i  : 0.0=dynamic ... 1.0=stable    (==st->env_stab )     */
                        Word16 *     alpha,           /*  o  : Magnitude modification factors for fade to average     */
                        Word16 *     beta,            /*    : Magnitude modification factors for fade to average     */
                        Word16 *     beta_mute, /* i/o  : Factor for long-term mute                              */
                        Word16 *     Xavg,      /* o  : Frequency group average gain to fade to                */
                     Word16 Q_spec,
                     Word32 L_oold_xfp_w_E_fx,
                     Word16 oold_xfp_w_E_exp_fx,
                     Word16 oold_Ltot_exp_fx,
                     Word16 *oold_grp_shape_fx,
                     Word32 L_old_xfp_w_E_fx,
                     Word16 old_xfp_w_E_exp_fx,
                     Word16 old_Ltot_exp_fx,
                     Word16 *old_grp_shape_fx,
                        Word8 *scratchBuffer);

void spec_ana_fx(Word16 *xfp, Word16 *, Word32 *, Word16 *, Word16 *, const Word16, const Word16,
                 const Word16 *, const Word16, const Word16, Word16 maxLprot, Word16 maxPlocs, Word8 *scratchBuffer);

void subst_spec_fx(const Word16 *, const Word32 *, Word16 *, const Word16, Word16 *, const Word16 *, const Word16,
                   const Word16 *, const Word16, Word16 *, const Word16 *, const Word16 *,
                   const Word16 *, const Word16
);

void rec_frame_fx(Word16 *     X,            /* i  : FFT spectrum */
                  Word32 *     L_ecu_rec,    /* o  : Reconstructed frame in tda domain */
                  const Word16 output_frame, /* i  : Frame length */
                  const Word16 Q,
                  const Word16 *const win2ms_init ,   /* i:  2 ms initial part of pre_tda window */
                  const Word16 *const win16ms_center, /* i:  16 ms combined part  of pre_tda IWHR+MDCT-ana  */

                  Word16 maxLprot,
                  const Word16 *prevsynth,
                  const Word16 Q_prevsynth,
                  Word8 *scratchBuffer);

Word16 rand_phase_fx(const Word16 seed, Word16 *sin_F, Word16 *cos_F);

void hq_phase_ecu_fx(const Word16 *prevsynth, /* i  : buffer of previously synthesized signal   */
                     Word32 *L_ecu_rec, /* o  : reconstructed frame in tda domain  , also tmp w32_fft buffer        */
                     Word16 *time_offs, /* i/o: Sample offset for consecutive frame losses*/
                     Word16 *X_sav,     /* i/o: Stored spectrum of prototype frame        */
                     Word16 *Q_spec,    /*  o: Q value of stored spectrum                */
                     Word16 *num_p,     /* i/o: Number of identified peaks                */
                     Word16 *plocs,     /* i/o: Peak locations                            */
                     Word32 *L_plocsi,  /* i/o: Interpolated peak locations           Q16 */
                     const Word16 env_stab,            /* i  : Envelope stability parameter              */
                     const Word16 f0hzLtpBinQ7,        /* i:  LTP bin frequency in normalized Hz  Q7 */
                     const Word16 norm_corrQ15_fx,     /*i : correlation for lag at f0hzLtpBinQ7 */
                     const Word16 prev_bfi,            /* i   : indicating burst frame error             */
                     Word16 old_is_transient[2], /* i/o   : flags indicating noise generation */
                     Word16 *     mag_chg_1st,         /* i/o: per band magnitude modifier for transients */
                     Word16 *     mag_chg_gr,          /*  o: per band magnitude modifier incl burst attenuation   */
                     Word16 *     Xavg,                /* i/o: Frequency group average gain to fade to   */
                     Word16 *     beta_mute,           /* o   : Factor for long-term mute                */
                     const Word16 bwidth_fx,           /* i  : Encoded bandwidth                         */
                     const Word16 output_frame,        /* i   : frame length                             */
                     Word16 * seed_out_fxPtr,           /* o: seed synch analysis */ 
                     Word16 * X_out ,                   /* o: utput  evolved spectrum  */  
                     const Word16  t_adv,      /* i  : time adjustment including time_offs       */
                     const Word16 *const win2ms_init ,   /* i:  2 ms initial part of pre_tda window */
                     const Word16 *const win16ms_center, /* i:  16 ms combined part  of pre_tda IWHR+MDCT-ana  */
                     const Word16 *sp_ana_win, /* i  : whr hamming window */
                     Word16 q_fx_old_exp, Word16 maxLprot,
                     Word16 maxPlocs,
                     Word32 L_oold_xfp_w_E_fx,
                     Word16 oold_xfp_w_E_exp_fx, /* exp of time signal */
                     Word16 oold_Ltot_exp_fx, /*true exp of energy */
                     Word16 *oold_grp_shape_fx,
                     Word32 L_old_xfp_w_E_fx,
                     Word16 old_xfp_w_E_exp_fx, /* exp of time signal */
                     Word16 old_Ltot_exp_fx,   /*true exp of energy */
                     Word16 *old_grp_shape_fx,
                     Word16  margin_prev_synth,   /* i: margin in prev_synth(16ms for first bfi , 3.75 ms for other bfi frames ) ,  from  plcAd.PhECU_margin_xfp */
                     Word8 *scratchBuffer /* Size = 2 * MAX_LGW + 8 * MAX_LPROT + 12 * MAX_L_FRAME */
);


Word16 plc_xcorr_lc_fx(/* o: quantized output xcorr in Q15  [ 0 ..32767 ] = [0. 1.0[  */
                Word16 *pcmbuf_fx, /* NB should be an  already dynamically upscaled buffer  with about 0...1  bits margin */
                Word16 buflen,    /* Q0 size of pcmbuf_fx */
                Word16  pitch_int, /* Q0  in Fs, lag value to evaluate, corresponding to the current f0    fr pcm_buf   */
                Word16 nom_corr_len, /* nominal correlation length to use for lags <= MAX_PITCH */
                Word16 fs_idx);
 
void plc_phEcu_peak_locator_fx(const Word16 *, const Word16, Word16 *, Word16 *, const Word16, const Word16, const Word16, Word16, Word8 *);

Word16  plc_phEcu_find_ind_fx(const Word16*  , const Word16   , const Word16   );


Word16 initQV(Word16 SR_idx, Word32 BR);

void processEstimateGlobalGain_fx(Word32 x[], Word16 x_e, Word16 lg, Word16 sqTargetBits, Word16 *gain, Word16 *gain_e,
                                  Word16 *quantizedGain, Word16 *quantizedGainMin, Word16 quantizedGainOff,
                                  Word32 *targetBitsOff, Word16 *old_targetBits, Word16 old_specBits,
                                  Word8 *scratchBuffer);

void processAdjustGlobalGain_fx(Word16 *gg_idx, Word16 gg_idx_min, Word16 gg_idx_off, Word16 *gain, Word16 *gain_e,
                                Word16 target, Word16 nBits, Word16 *gainChange, Word16 fs_idx);

void processScalarQuant_fx(Word32 x[], Word16 x_e, Word16 xq[], Word16 L_frame, Word16 gain, Word16 gain_e);

void processQuantizeSpec_fx(Word32 x[], Word16 x_e, Word16 gain, Word16 gain_e, Word16 xq[], Word16 nt, Word16 target,
                            Word16 totalBits, Word16 *nBits, Word16 *nBits2, Word16 fs_idx, Word16 *lastnz,
                            Word16 *codingdata, Word16 *lsbMode, Word16 mode);

void processResidualCoding_fx(Word16 x_e, Word32 x[], Word16 xq[], Word16 gain, Word16 gain_e, Word16 L_spec,
                              Word16 targetBits, Word16 nBits, Word8 *resBits, Word16 *numResBits);

void scale_signal24_fx(Word32 const x[], /* i:   time input signal */
                       Word16 x_scaled[], Word16 *x_exp, Word16 mdct_mem[], Word16 mdct_mem_len,
                       Word16 resample_mem_in[], Word16 resample_mem_in_len, Word32 resample_mem_in50[],
                       Word16 resample_mem_out[], Word16 resample_mem_out_len, Word32 mdct_mem32[], Word16 N,
                       Word32 resamp_mem32[], Word16 mem_s12k8[], Word16 *resamp_scale);

void processReorderBitstream_fx(UWord8 *bytes, Word16 n_pccw, Word16 n_pc, Word16 b_left, Word8 *scratchBuffer);

/* al_fec.c */
Word16 fec_get_n_pccw(Word16 slot_bytes, Word16 fec_mode, Word16 ccc_flag);
Word16 fec_get_data_size(Word16 fec_mode, Word16 ccc_flag, Word16 slot_bytes);
Word16 fec_get_n_pc(Word16 fec_mode, Word16 n_pccw, Word16 slot_bytes);

void fec_encoder(Word16 mode, Word16 epmr, UWord8 *iobuf, Word16 data_bytes, Word16 slot_bytes, Word16 n_pccw,
                 void *scratch);

int fec_decoder(UWord8 *iobuf, Word16 slot_bytes, int *data_bytes, Word16 *epmr, Word16 ccc_flag, Word16 *n_pccw,
                int *bfi, Word16 *be_bp_left, Word16 *be_bp_right, Word16 *n_pc, Word16 *m_fec, void *scratch);

void processPCmain_fx(Word16 rframe, Word16 *bfi, Word16 prev_bfi, Word16 yLen, Word16 frame_dms, Word16 q_old_res_fx[],
                      Word16 *q_old_res_fx_exp, Word16 q_res_fx[], Word16 q_old_d_fx[], Word16 spec_inv_idx,
                      Word16 pitch_present, Word16 stab_fac, Word32 q_d_fx[], Word16 *q_fx_exp,
                      Word16 gg_idx, Word16 gg_idx_off, Word16 *prev_gg, Word16 *prev_gg_e, Word16 *BW_cutoff_idx_nf,
                      Word16 *prev_BW_cutoff_idx_nf, Word16 fac_ns_idx, Word16 *prev_fac_ns_fx, Word16 *pc_nbLostFramesInRow);
void processPCclassify_fx(Word16 pitch_present, Word16 frame_dms, Word16 q_old_d_fx[], Word16 q_old_res_fx[],
                          Word16 yLen, Word16 spec_inv_idx, Word16 stab_fac, Word16 prev_bfi, Word16 *bfi);
void processPCapply_fx(Word16 yLen, Word16 q_old_res_fx[], Word16 q_res_fx[], Word16 *q_old_res_fx_exp, Word16 q_old_d_fx[],
                       Word16 spec_inv_idx, Word16 *fac, Word16 *fac_e, Word32 q_d_fx[], Word16 *q_fx_exp,
                       Word16 gg_idx, Word16 gg_idx_off, Word16 prev_gg, Word16 prev_gg_e, Word16 *pc_nbLostFramesInRow);
void processPCupdate_fx(Word16 bfi, Word16 yLen, Word16 q_old_res_fx[], Word16 *q_old_res_fx_exp,
                        Word16 q_res_fx[], Word16 spec_inv_idx, Word16 gg_idx, Word16 gg_idx_off,
                        Word16 *prev_gg, Word16 *prev_gg_e, Word16 rframe, Word16 *BW_cutoff_idx_nf,
                        Word16 *prev_BW_cutoff_idx_nf, Word16 fac_ns_idx, Word16 *prev_fac_ns_fx, Word16 fac, Word16 fac_e);
void processPcApplyDamping_fx(Word32 x[], Word16 xLen, Word16 fac, Word16 spec_inv_idx);

void process_cutoff_bandwidth(Word32 d_fx[], Word16 len, Word16 bw_bin);

void dct16_fx(const Word16 *in, Word16 *out);
void idct16_fx(const Word16 *in, Word16 *out);

/* Functions used in arithmetic coder */

void write_bit_backward(UWord8 *ptr, Word16 *bp, Word16 *mask, Word16 bit);
void write_indice_backward(UWord8 *ptr, Word16 *bp, Word16 *mask, Word16 indice, Word16 numbits);

void processEncoderEntropy(UWord8 *bytes, Word16 *bp_side, Word16 *mask_side, Word16 nbbits, Word16 targetBytes,
                           Word16 L_spec, Word16 BW_cutoff_bits, Word16 tns_numfilters,
                           Word16 lsbMode, Word16 lastnz, Word16 *tns_order, Word16 fac_ns_idx, Word16 gg_idx,
                           Word16 BW_cutoff_idx, Word16 *ltpf_idx, Word32 *L_scf_idx, Word16 bfi_ext, Word16 fs_idx);

void processDecoderEntropy_fx(UWord8 *bytes, Word16 *bp_side, Word16 *mask_side, Word16 nbbits,
                              Word16 L_spec, Word16 fs_idx, Word16 BW_cutoff_bits, Word16 *tns_numfilters,
                              Word16 *lsbMode, Word16 *lastnz, Word16 *bfi, Word16 *tns_order, Word16 *fac_ns_idx,
                              Word16 *gg_idx, Word16 *BW_cutoff_idx, Word16 *ltpf_idx, Word32 *L_scf_idx,
                              Word16 frame_dms);

#ifdef ENABLE_PADDING
int paddingDec_fx(UWord8 *bytes, Word16 nbbits, Word16 L_spec, Word16 BW_cutoff_bits, Word16 ep_enabled,
                   Word16 *total_padding, Word16 *np_zero);
#endif

Word16 read_bit(UWord8 *ptr, Word16 *bp, Word16 *mask);

/* setup_enc_lc3.c */
int       alloc_encoder(LC3_Enc *encoder, int samplerate, int channels);
void      set_enc_frame_params(LC3_Enc *encoder);
LC3_Error update_enc_bitrate(LC3_Enc *encoder, int bitrate);
LC3_Error FillEncSetup(LC3_Enc *encoder, int samplerate, int channels);

/* setup_dec_lc3.c */
int       alloc_decoder(LC3_Dec *decoder, int samplerate, int channels, LC3_PlcMode plc_mode);
void      set_dec_frame_params(LC3_Dec *decoder);
LC3_Error update_dec_bitrate(LC3_Dec *decoder, int ch, Word16 nBytes);
LC3_Error FillDecSetup(LC3_Dec *decoder, int samplerate, int channels, LC3_PlcMode plc_mode);

int       Enc_LC3(LC3_Enc *encoder, void **input, int bits_per_sample, UWord8 *output, void *scratch, Word16 bfi_ext);
LC3_Error Dec_LC3(LC3_Dec *decoder, UWord8 *input, int input_bytes, void **output, int bits_per_sample, void *scratch, int bfi_ext);

void *balloc(void *base, size_t *base_size, size_t size);


#endif
