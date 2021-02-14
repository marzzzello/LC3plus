/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef SETUP_DEC_LC3_H
#define SETUP_DEC_LC3_H

#include "constants.h"

typedef struct
{
    Word16 *x_old_tot_fx;      /* MAX_LEN_PCM_PLC_TOT    */
    Word32 *PhECU_f0est;       /* MAX_PLOCS            interpolated plocs  */
    Word16 *PhECU_xfp_fx;      /* MAX_LPROT */
    Word16 *PhECU_X_sav_fx;    /* MAX_LPROT */
    Word16 *PhECU_plocs;       /* MAX_PLOCS */
    Word16 *PhECU_fg_wintaper; /* MDCT_MEM_LEN_MAX */
    Word16 *PhECU_win_pre_tda; /* MAX_WIN_PRE_TDA */
    Word16  tdc_gain_p;
    Word32  tdc_gain_c;
    Word16  stab_fac;
    Word16  tdc_fract;
    Word16  tdc_seed;
    Word16  tdc_cum_damp;
    Word16  tdc_preemph_fac;
    Word16  tdc_lpc_order;
    Word16  cum_fading_slow;
    Word16  cum_fading_fast;
    Word16  PhECU_LprotOrg_fx; /* needed to change the Prot size  adaptively  */
    Word16  PhECU_Lprot_fx;
    Word16  PhECU_fs_idx_fx;
    Word16  PhECU_frame_ms;  /* needed in PLC_Update and PLCMain functons*/  
    Word16  PhECU_seed_fx;
    Word16  PhECU_xfp_exp_fx;
    Word16  PhECU_time_offs;
    Word16  PhECU_X_savQ_fx;
    Word16  PhECU_num_plocs;
    Word16  PhECU_f0hzLtpBinQ7; /*  ltp F0 in bins  if available  */
    Word16  PhECU_short_flag_prev;
    Word16  PhECU_whr_tot_taper;
    Word16  PhECU_whr_tot_flat;
    Word16  PhECU_LDWIN_OLAP;
    Word16  PhECU_LA;
    Word16  PhECU_t_adv;
    Word16  PhECU_beta_mute;
    Word16  norm_corrQ15_fx;
    Word16  q_fx_old_exp;
    Word16  max_len_pcm_plc;
    Word16  max_lprot;
    Word16  max_plocs;
    
    /* Word32 L_tot W_energy sum exponent */ 
    Word16  PhECU_oold_Ltot_exp_fx; 
    Word16  PhECU_old_Ltot_exp_fx;
    Word32  PhECU_L_oold_xfp_w_E_fx;
    Word32  PhECU_L_old_xfp_w_E_fx;
    Word16  PhECU_oold_xfp_w_E_exp_fx;   /* input Word16 xfp exponnet  */
    Word16  PhECU_old_xfp_w_E_exp_fx;  
    Word16  PhECU_oold_grp_shape_fx[MAX_LGW] ALIGN_BUFFER_STRUCT;
    Word16  PhECU_old_grp_shape_fx[MAX_LGW] ALIGN_BUFFER_STRUCT;
    Word16  PhECU_margin_xfp; 
    Word16  PhECU_mag_chg_1st[MAX_LGW] ALIGN_BUFFER_STRUCT;
    Word16  PhECU_Xavg[MAX_LGW] ALIGN_BUFFER_STRUCT;
    Word16  old_scf_q[M] ALIGN_BUFFER_STRUCT;
    Word16  old_old_scf_q[M] ALIGN_BUFFER_STRUCT;
    Word16  tdc_A[M + 1] ALIGN_BUFFER_STRUCT;
    /* for now 20 ms saved Q14  or ptr to a combined ifft win and MDCT  preTDA synthesis window  16  ms */
} AplcSetup;

/* Channel state and bitrate-derived values go in this struct */
typedef struct
{
    Word16 *ltpf_mem_x;       /* LTPF_MEM_X_LEN */
    Word16 *ltpf_mem_y;       /* LTPF_MEM_Y_LEN */
    Word16 *stDec_ola_mem_fx; /* MDCT_MEM_LEN_MAX */
    AplcSetup *plcAd;
    Word16 *   q_old_d_fx; /* MAX_BW */
    Word16     q_old_fx_exp;
    Word16     ns_seed;
    Word16     ns_cum_alpha;
    Word16  pc_seed;
    Word16  pc_nbLostFramesInRow;
    Word16 *q_old_res_fx;
    Word16  q_old_res_fx_exp;
    Word16  prev_gg;
    Word16  prev_gg_e;
    Word16  prev_BW_cutoff_idx_nf;
    Word16 prev_fac_ns_fx;
    Word16 total_bits;
    Word16 enable_lpc_weighting;
    Word16 stDec_ola_mem_fx_exp;
    Word16 targetBytes;
    Word16 ltpf_mem_e;
    Word16 ltpf_mem_pitch_int;
    Word16 ltpf_mem_pitch_fr;
    Word16 ltpf_mem_gain;
    Word16 ltpf_mem_active;
    Word16 ltpf_scale_fac_idx;
    Word16 ltpf_mem_scale_fac_idx;
    Word16 quantizedGainOff;
    Word16 prev_bfi;
    Word16 prev_prev_bfi;
    Word16 concealMethod;
    Word16 nbLostFramesInRow;
    Word16 plc_damping;
    Word16 last_size;
} DecSetup;

/* Constants and sampling rate derived values go in this struct */
struct LC3_Dec
{
    DecSetup *    channel_setup[MAX_CHANNELS];
    const Word16 *W_fx;
    const Word16 *bands_offset;
    Word32        fs;           /* sampling rate, 44.1 maps to 48 */
    Word32        fs_out;       /* output sampling rate */
    Word16        fs_idx;       /* sampling rate index */
    Word16        frame_length; /* sampling rate index */
    Word16        channels;     /* number of channels */
    Word16        plcMeth;      /* PLC method for all channels */
    Word16        frame_dms;    /* frame length in dms (decimilliseconds, 10^-4)*/
    Word16        last_size;    /* size of last frame, without error protection */
    Word16        ep_enabled;   /* error protection enabled */
    Word16        error_report; /* corrected errors in last frame or -1 on error */

    Word16 n_pccw;
    Word16 be_bp_left;
    Word16 be_bp_right;
    Word16 n_pc;
    Word16 m_fec;
    Word16 epmr;
    Word16 combined_channel_coding;

    Word16 yLen;
    Word16 W_size;
    Word16 la_zeroes;
    Word16 stDec_ola_mem_fx_len;
    Word16 bands_number;
    Word16 ltpf_mem_x_len;
    Word16 ltpf_mem_y_len;
    Word16 BW_cutoff_bits;
};

#endif
