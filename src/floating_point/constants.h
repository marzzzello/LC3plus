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
#include "structs.h"

/* Ari coder */
extern const LC3_INT ari_tns_order_cf[2][9];
extern const LC3_INT ari_tns_freq_cf[8][18];
extern const LC3_INT ari_spec_lookup_fl[4096];
extern const LC3_INT ari_spec_cumfreq_fl[64][18];
extern const LC3_INT ari_spec_bits_fl[64][17];

/* SNS */
extern const LC3_FLOAT sns_C1[8][32];
extern const LC3_FLOAT sns_C2[8][32];
extern const LC3_INT   pvq_enc_A[16][11];

/* 12.8 kHz resampler */
extern const LC3_FLOAT lp_scale_factors[6];
extern const LC3_FLOAT lp_filter[240];
extern const double highpass50_filt_b[3];
extern const double highpass50_filt_a[3];
extern const LC3_INT    up_fac[6];

/* TNS */
extern const LC3_FLOAT quants_pts_tns[17];
extern const LC3_INT   huff_bits_tns[8][17];
extern const LC3_INT   order1_tns[8];
extern const LC3_INT   order2_tns[8];
extern const LC3_FLOAT lagw_tns[9];
extern const LC3_FLOAT quants_pts_tns[17];
extern const LC3_FLOAT quants_thr_tns[18];

/* SNS */
extern const LC3_FLOAT sns_vq_far_adj_gains_fl[8];
extern const LC3_FLOAT sns_vq_near_adj_gains_fl[4];
extern const LC3_FLOAT sns_vq_reg_lf_adj_gains_fl[4];
extern const LC3_FLOAT q_g_sns[6];
extern const LC3_FLOAT sns_vq_reg_adj_gains_fl[2];
extern const LC3_FLOAT sns_dec_gains[4][8];

/* Global Gain */
extern const LC3_INT   gg_p1[6];
extern const LC3_INT   gg_p2[6];
extern const LC3_INT   gg_p3[6];
extern const LC3_FLOAT gg_c[6];
extern const LC3_FLOAT gg_d[6];

/* Olpa */
extern const LC3_FLOAT olpa_down2[5];
extern const LC3_FLOAT olpa_acw[98];

/* LTPF */
extern const LC3_FLOAT conf_inter_filter_48[4][12];
extern const LC3_FLOAT conf_inter_filter_32[4][8];
extern const LC3_FLOAT conf_inter_filter_24[4][6];
extern const LC3_FLOAT conf_inter_filter_16[4][4];
extern const LC3_FLOAT conf_tilt_filter_48[4][11];
extern const LC3_FLOAT conf_tilt_filter_32[4][7];
extern const LC3_FLOAT conf_tilt_filter_24[4][5];
extern const LC3_FLOAT conf_tilt_filter_16[4][3];
extern const LC3_FLOAT inter4_1[33];
extern const LC3_FLOAT enc_inter_filter[4][4];

/* Bandwidth Detector */
extern const LC3_INT  threshold_quiet[4];
extern const LC3_INT  threshold_brickwall[4];
extern const LC3_INT  brickwall_dist[4];
extern const LC3_INT  BW_warp_idx_start_16k[4];
extern const LC3_INT  BW_warp_idx_stop_16k[4];
extern const LC3_INT  BW_warp_idx_start_24k[4];
extern const LC3_INT  BW_warp_idx_stop_24k[4];
extern const LC3_INT  BW_warp_idx_start_32k[4];
extern const LC3_INT  BW_warp_idx_stop_32k[4];
extern const LC3_INT  BW_warp_idx_start_48k[4];
extern const LC3_INT  BW_warp_idx_stop_48k[4];
extern const LC3_INT* BW_warp_idx_start_all[4];
extern const LC3_INT* BW_warp_idx_stop_all[4];

extern const LC3_INT  BW_warp_idx_start_16k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_16k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_start_24k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_24k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_start_32k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_32k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_start_48k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_48k_2_5ms[4];
extern const LC3_INT* BW_warp_idx_start_all_2_5ms[4];
extern const LC3_INT* BW_warp_idx_stop_all_2_5ms[4];
#ifdef ENABLE_HR_MODE
extern const LC3_INT  BW_cutoff_bin_all_2_5ms_HR[MAX_BW_BANDS_NUMBER];
extern const LC3_INT  BW_cutoff_bin_all_2_5ms[MAX_BW_BANDS_NUMBER];
extern const LC3_INT  bands_number_2_5ms_HR[6];
extern const LC3_INT  bands_number_2_5ms[6];
#else
extern const LC3_INT  BW_cutoff_bin_all_2_5ms[MAX_BW_BANDS_NUMBER];
extern const LC3_INT  bands_number_2_5ms[5];
#endif


extern const LC3_INT  BW_warp_idx_start_16k_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_16k_5ms[4];
extern const LC3_INT  BW_warp_idx_start_24k_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_24k_5ms[4];
extern const LC3_INT  BW_warp_idx_start_32k_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_32k_5ms[4];
extern const LC3_INT  BW_warp_idx_start_48k_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_48k_5ms[4];
extern const LC3_INT* BW_warp_idx_start_all_5ms[4];
extern const LC3_INT* BW_warp_idx_stop_all_5ms[4];
extern const LC3_INT  BW_cutoff_bin_all_5ms[MAX_BW_BANDS_NUMBER];
#ifdef ENABLE_HR_MODE
extern const LC3_INT bands_number_5ms[6];
extern const LC3_INT BW_cutoff_bin_all_HR[MAX_BW_BANDS_NUMBER];
extern const LC3_INT BW_cutoff_bin_all_5ms_HR[MAX_BW_BANDS_NUMBER];
#else
extern const LC3_INT    bands_number_5ms[5];
#endif
extern const LC3_INT BW_cutoff_bin_all[MAX_BW_BANDS_NUMBER];
extern const LC3_INT BW_cutoff_bits_all[MAX_BW_BANDS_NUMBER];

/* Arithmetic coding */
extern const LC3_INT tns_cf[8][18];
extern const LC3_INT tns_freq_cf[2][9];

/* MDCT Windows */
extern const LC3_FLOAT MDCT_WINDOW_80[160];
extern const LC3_FLOAT MDCT_WINDOW_160[320];
extern const LC3_FLOAT MDCT_WINDOW_240[480];
extern const LC3_FLOAT MDCT_WINDOW_320[640];
extern const LC3_FLOAT MDCT_WINDOW_480[960];
extern const LC3_FLOAT MDCT_WINDOW_960[1920];
#ifdef ENABLE_HR_MODE
extern const LC3_FLOAT* MDCT_WINS_10ms[2][6];
#else
extern const LC3_FLOAT* MDCT_WINS_10ms[5];
#endif

extern const LC3_FLOAT MDCT_WINDOW_80_2_5ms[40];
extern const LC3_FLOAT MDCT_WINDOW_160_2_5ms[80];
extern const LC3_FLOAT MDCT_WINDOW_240_2_5ms[120];
extern const LC3_FLOAT MDCT_WINDOW_320_2_5ms[160];
extern const LC3_FLOAT MDCT_WINDOW_480_2_5ms[240];
#ifdef ENABLE_HR_MODE
extern const LC3_FLOAT* MDCT_WINS_2_5ms[2][6];
#else
extern const LC3_FLOAT* MDCT_WINS_2_5ms[5];
#endif


extern const LC3_FLOAT MDCT_WINDOW_80_5ms[80];
extern const LC3_FLOAT MDCT_WINDOW_160_5ms[160];
extern const LC3_FLOAT MDCT_WINDOW_240_5ms[240];
extern const LC3_FLOAT MDCT_WINDOW_320_5ms[320];
extern const LC3_FLOAT MDCT_WINDOW_480_5ms[480];
#ifdef ENABLE_HR_MODE
extern const LC3_FLOAT* MDCT_WINS_5ms[2][6];
#else
extern const LC3_FLOAT* MDCT_WINS_5ms[5];
#endif

#ifdef ENABLE_HR_MODE
extern const LC3_INT MDCT_WINDOWS_LENGTHS_10ms[6];
#else
extern const LC3_INT    MDCT_WINDOWS_LENGTHS_10ms[5];
#endif

#ifdef ENABLE_HR_MODE
extern const LC3_INT    MDCT_WINDOWS_LENGTHS_2_5ms[6];
#else
extern const LC3_INT    MDCT_WINDOWS_LENGTHS_2_5ms[5];
#endif


#ifdef ENABLE_HR_MODE
extern const LC3_INT MDCT_WINDOWS_LENGTHS_5ms[6];
#else
extern const LC3_INT    MDCT_WINDOWS_LENGTHS_5ms[5];
#endif

/* Per band energy */
#ifdef ENABLE_HR_MODE
extern const LC3_INT* ACC_COEFF_PER_BAND[6];
extern const LC3_INT* ACC_COEFF_PER_BAND_HR[6];
#else
extern const LC3_INT*   ACC_COEFF_PER_BAND[6];
#endif
#ifdef ENABLE_HR_MODE
extern const LC3_INT* ACC_COEFF_PER_BAND_2_5ms[6];
extern const LC3_INT* ACC_COEFF_PER_BAND_2_5ms_HR[6];
#else
extern const LC3_INT*   ACC_COEFF_PER_BAND_2_5ms[5];
#endif
#ifdef ENABLE_HR_MODE
extern const LC3_INT* ACC_COEFF_PER_BAND_5ms[6];
extern const LC3_INT* ACC_COEFF_PER_BAND_5ms_HR[6];
#else
extern const LC3_INT*   ACC_COEFF_PER_BAND_5ms[5];
#endif



#endif
