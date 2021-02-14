/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "defines.h"
#include "basop_util.h"


extern RAM_ALIGN const Word16 LowDelayShapes_n960_len[5];
extern RAM_ALIGN const Word16 LowDelayShapes_n960_la_zeroes[5];
extern RAM_ALIGN const Word16 *const LowDelayShapes_n960[5];
extern RAM_ALIGN const Word16 LowDelayShapes_n960_len_5ms[5];
extern RAM_ALIGN const Word16 LowDelayShapes_n960_la_zeroes_5ms[5];
extern RAM_ALIGN const Word16 *const LowDelayShapes_n960_5ms[5];
extern RAM_ALIGN const Word16 LowDelayShapes_n960_len_2_5ms[5];
extern RAM_ALIGN const Word16 LowDelayShapes_n960_la_zeroes_2_5ms[5];
extern RAM_ALIGN const Word16 *const LowDelayShapes_n960_2_5ms[5];

extern RAM_ALIGN const Word32 BW_thresh_quiet[4];
extern RAM_ALIGN const Word16 BW_thresh_quiet_exp;
extern RAM_ALIGN const Word16 BW_thresh_brickwall[4];
extern RAM_ALIGN const Word16 BW_brickwall_dist[4];
extern RAM_ALIGN const Word16 BW_cutoff_bin_all[MAX_BW_BANDS_NUMBER];
extern RAM_ALIGN const Word16 BW_cutoff_bits_all[MAX_BW_BANDS_NUMBER];
extern RAM_ALIGN const Word16 *const BW_warp_idx_start_all[MAX_BW_BANDS_NUMBER - 1];
extern RAM_ALIGN const Word16 *const BW_warp_idx_stop_all[MAX_BW_BANDS_NUMBER - 1];
extern RAM_ALIGN const Word16 *const BW_warp_idx_start_all_5ms[MAX_BW_BANDS_NUMBER - 1];
extern RAM_ALIGN const Word16 *const BW_warp_idx_stop_all_5ms[MAX_BW_BANDS_NUMBER - 1];
extern RAM_ALIGN const Word16 *const BW_warp_idx_start_all_2_5ms[MAX_BW_BANDS_NUMBER - 1];
extern RAM_ALIGN const Word16 *const BW_warp_idx_stop_all_2_5ms[MAX_BW_BANDS_NUMBER - 1];

extern RAM_ALIGN const Word16 *const tns_subdiv_startfreq[MAX_BW_BANDS_NUMBER];
extern RAM_ALIGN const Word16 *const tns_subdiv_stopfreq[MAX_BW_BANDS_NUMBER];
extern RAM_ALIGN const Word16 *const tns_subdiv_startfreq_5ms[MAX_BW_BANDS_NUMBER];
extern RAM_ALIGN const Word16 *const tns_subdiv_stopfreq_5ms[MAX_BW_BANDS_NUMBER];
extern RAM_ALIGN const Word16 *const tns_subdiv_startfreq_2_5ms[MAX_BW_BANDS_NUMBER];
extern RAM_ALIGN const Word16 *const tns_subdiv_stopfreq_2_5ms[MAX_BW_BANDS_NUMBER];
extern RAM_ALIGN const Word16 Tab_esc_nb[4];

extern RAM_ALIGN const Word8 ari_spec_lookup[4096];
extern RAM_ALIGN const UWord16 ari_spec_cumfreq[64][17];
extern RAM_ALIGN const UWord16 ari_spec_freq[64][17];
extern RAM_ALIGN const UWord16 ari_spec_bits[64][17];

extern RAM_ALIGN const Word32 tnsAcfWindow[MAXLAG];
extern RAM_ALIGN const Word16 ac_tns_order_bits[2][MAXLAG];
extern RAM_ALIGN const Word16 ac_tns_order_freq[2][MAXLAG];
extern RAM_ALIGN const Word16 ac_tns_order_cumfreq[2][MAXLAG];
extern RAM_ALIGN const Word16 ac_tns_coef_bits[MAXLAG][TNS_COEF_RES];
extern RAM_ALIGN const Word16 ac_tns_coef_freq[MAXLAG][TNS_COEF_RES];
extern RAM_ALIGN const Word16 ac_tns_coef_cumfreq[MAXLAG][TNS_COEF_RES];
extern RAM_ALIGN const Word16 tnsQuantPts[TNS_COEF_RES];
extern RAM_ALIGN const Word16 tnsQuantThr[TNS_COEF_RES - 1];

extern RAM_ALIGN const Word16 *const lpc_pre_emphasis[5];
extern RAM_ALIGN const Word16 *const lpc_pre_emphasis_e[5];

extern RAM_ALIGN const Word16 *const lpc_lin_pre_emphasis[5];
extern RAM_ALIGN const Word16 *const lpc_lin_pre_emphasis_e[5];
extern RAM_ALIGN const Word16 *const lpc_lin_pre_emphasis_5ms[5];
extern RAM_ALIGN const Word16 *const lpc_lin_pre_emphasis_e_5ms[5];
extern RAM_ALIGN const Word16 *const lpc_lin_pre_emphasis_2_5ms[5];
extern RAM_ALIGN const Word16 *const lpc_lin_pre_emphasis_e_2_5ms[5];

extern RAM_ALIGN const Word16 *const lpc_warp_dee_emphasis[5];
extern RAM_ALIGN const Word16 *const lpc_warp_dee_emphasis_e[5];

extern RAM_ALIGN const Word16 bands_nrg_scale[32];

extern RAM_ALIGN const Word16 *const bands_offset[5];
extern RAM_ALIGN const Word16 bands_offset_with_one_max[5];
extern RAM_ALIGN const Word16 bands_offset_with_two_max[5];
extern RAM_ALIGN const Word16 bands_number_5ms[5];
extern RAM_ALIGN const Word16 *const bands_offset_5ms[5];
extern RAM_ALIGN const Word16 bands_offset_with_one_max_5ms[5];
extern RAM_ALIGN const Word16 bands_offset_with_two_max_5ms[5];
extern RAM_ALIGN const Word16 bands_number_2_5ms[5];
extern RAM_ALIGN const Word16 *const bands_offset_2_5ms[5];
extern RAM_ALIGN const Word16 bands_offset_with_one_max_2_5ms[5];
extern RAM_ALIGN const Word16 bands_offset_with_two_max_2_5ms[5];

extern RAM_ALIGN const Word16 pitch_max[5];
extern RAM_ALIGN const Word16 plc_preemph_fac[5];

extern RAM_ALIGN const Word16 TDC_high_16[11];
extern RAM_ALIGN const Word16 TDC_high_32[11];

extern RAM_ALIGN const Word32 *const lag_win[5];

extern RAM_ALIGN const Word16 *const bands_offset_lin[5];
extern RAM_ALIGN const Word16 bands_offset_with_one_max_lin[5];
extern RAM_ALIGN const Word16 bands_offset_with_two_max_lin[5];
extern RAM_ALIGN const Word16 *const bands_offset_lin_5ms[5];
extern RAM_ALIGN const Word16 bands_offset_with_one_max_lin_5ms[5];
extern RAM_ALIGN const Word16 bands_offset_with_two_max_lin_5ms[5];
extern RAM_ALIGN const Word16 *const bands_offset_lin_2_5ms[5];
extern RAM_ALIGN const Word16 bands_offset_with_one_max_lin_2_5ms[5];
extern RAM_ALIGN const Word16 bands_offset_with_two_max_lin_2_5ms[5];

extern RAM_ALIGN const Word32 inv_odft_twiddle_80_re[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_80_im[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_60_re[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_60_im[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_40_re[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_40_im[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_20_re[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_20_im[M];

#ifdef SUBSET_SQ
extern RAM_ALIGN const Word16 resamp_filt_16k[240];
#else
extern RAM_ALIGN const Word16 resamp_filt_16k[1];
#endif
#ifdef SUBSET_HQ
extern RAM_ALIGN const Word16 resamp_filt_24k[240];
#else
extern RAM_ALIGN const Word16 resamp_filt_24k[1];
#endif
#ifdef SUBSET_SWB
extern RAM_ALIGN const Word16 resamp_filt_32k[240];
#else
extern RAM_ALIGN const Word16 resamp_filt_32k[1];
#endif
#ifdef SUBSET_FB
extern RAM_ALIGN const Word16 resamp_filt_48k[240];
#else
extern RAM_ALIGN const Word16 resamp_filt_48k[1];
#endif

extern RAM_ALIGN const Word16 resamp_params[5][4];
extern RAM_ALIGN const Word16 *const resamp_filts[5];

extern RAM_ALIGN const Word16 highpass50_filt_num[3];
extern RAM_ALIGN const Word16 highpass50_filt_den[2];

extern RAM_ALIGN const Word16 olpa_ac_weighting[98];

extern RAM_ALIGN const Word16 ltpf_ac_interp_filt[7][9];
extern RAM_ALIGN const Word16 inter_filter[5][4][12];
extern RAM_ALIGN const Word16 inter_filter_shift[5];
extern RAM_ALIGN const Word16 inter_filter_len[5];
extern RAM_ALIGN const Word16 tilt_filter[5][4][11];
extern RAM_ALIGN const Word16 tilt_filter_len[5];
extern RAM_ALIGN const Word16 gain_scale_fac[4];
extern RAM_ALIGN const Word16 pitch_scale[5];

typedef struct
{
    Word16  lead_sign_ind; /* this MPVQ index  is the first  part   of the total PVQ index  */
    UWord32 index, size;   /* this MPVQ index  is the second part  of the total PVQ index */
    Word16  dim, k_val;
    Word16  vec[PVQ_MAX_VEC_SIZE]; /* integer vector */
} PvqEntry_fx;

extern RAM_ALIGN const Word16 sns_vq_reg_adj_scf[2];
extern RAM_ALIGN const Word16 sns_vq_reg_lf_adj_scf[4];
extern RAM_ALIGN const Word16 sns_vq_near_adj_scf[4];
extern RAM_ALIGN const Word16 sns_vq_far_adj_scf[8];
extern RAM_ALIGN const Word16 *const sns_gaintabPtr[4];
extern RAM_ALIGN const Word16 sns_gainSz[4];
extern RAM_ALIGN const Word16 sns_gainMSBbits[4];
extern RAM_ALIGN const Word16 sns_gainLSBbits[4];
extern RAM_ALIGN const Word16 sns_Kval[4][2];
extern RAM_ALIGN const UWord32 sns_MPVQ_Sz[4][2];

extern RAM_ALIGN const Word16 st1SCF0_7_base5_32x8_Q14[256];
extern RAM_ALIGN const Word16 st1SCF8_15_base5_32x8_Q14[256];

/* PVQ deindexing tables */
extern RAM_ALIGN const UWord32 h_memN16K12[12 + 2];
extern RAM_ALIGN const UWord32 h_memN10K22[22 + 2];
extern RAM_ALIGN const UWord32 h_memN6K2[2 + 2];
extern RAM_ALIGN const Word16 tabledKMAX[16 + 1];
extern RAM_ALIGN const UWord32 *const MPVQ_offs_ptr[16 + 1];

extern RAM_ALIGN const Word16 isqrt_Q16tab[1 + SQRT_EN_MAX_FX];

extern RAM_ALIGN const Word16 adjust_global_gain_tables[5][5];


extern RAM_ALIGN const Word16 sqrt_table_phecu[];
extern RAM_ALIGN const Word16  POW_ATT_TABLE0[];
extern RAM_ALIGN const Word16  POW_ATT_TABLE1[];
#ifdef PLC2_FADEOUT_IN_MS 
#if PLC2_FADEOUT_IN_MS == 0
extern RAM_ALIGN const Word16* const POW_ATT_TABLES[3];
#else
extern RAM_ALIGN const Word16* const POW_ATT_TABLES[11];
#endif
#else
extern RAM_ALIGN const Word16* const POW_ATT_TABLES[3];
#endif

extern RAM_ALIGN const Word16 e_tot_headroom[];
extern RAM_ALIGN const Word16 xfp_wE_MDCT2FFTQ11[];
 
extern RAM_ALIGN const Word16 num_FsByResQ0[5];
extern RAM_ALIGN const Word16* const LprotSzPtr;  
extern RAM_ALIGN const   Word16  InvLprot_Q22[5];
extern RAM_ALIGN const   Word16 PhEcuFftScale[5];
#ifdef  NONBE_PLC2_MUTING_DCSYNT_FIX 
extern RAM_ALIGN const    Word16  oneOverFrameQ15Tab[5];
#endif
extern RAM_ALIGN const    Word16 PhEcu_Xsav_Flt2FxDnShift[];
extern RAM_ALIGN const    Word16 PhEcu_Xsav_Flt2FxScaleQ15[];  
extern RAM_ALIGN const    Word16 PhEcu_frac_thr_rise_lin_Q15[];
extern RAM_ALIGN const    Word16 PhEcu_frac_thr_decay_lin_Q15[];

extern RAM_ALIGN const Word16 mdct_grp_bins_fx[];
extern RAM_ALIGN const Word16 xavg_N_grp_fx[];
extern RAM_ALIGN const Word16 spec_shape_headroom[];
extern RAM_ALIGN const   Word16 rectLengthTab[];
extern RAM_ALIGN const   Word16 hamm_len2Tab[];

extern RAM_ALIGN const Word16 gw_len_inv_shift_fx[];
extern RAM_ALIGN const Word16 gwlpr_fx[];

extern RAM_ALIGN const Word16 sin_quarterQ15_fx[];
extern RAM_ALIGN const Word16 sincos_lowres_tab_sinQ15_fx[];

extern RAM_ALIGN const Word16 *const PhECU_wins[5][3];

extern   RAM_ALIGN const  Word16 *const w_new[];
extern   RAM_ALIGN const  Word16 *const w_old[];

/* extern    RAM_ALIGN const  Word16 WORK_LEN[]; */
extern    RAM_ALIGN const  Word16 COPY_LEN[];
extern    RAM_ALIGN const Word16 OLA_LEN[];


#endif
