/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef SETUP_DEC_LC3_FL_H
#define SETUP_DEC_LC3_FL_H

#include "constants.h"

/* Channel state and bitrate-derived values go in this struct */
typedef struct {
    LC3_INT* stDec_ola_mem_fx; /* MDCT_MEM_LEN_MAX */
    LC3_INT  total_bits;
    LC3_INT  enable_lpc_weighting;
    LC3_INT  targetBytes;
    LC3_INT  ltpf_mem_active;
    LC3_INT  quantizedGainOff;
    LC3_INT  nbLostFramesInRow;
    LC3_INT  ltpf_param[3];
    LC3_INT  ltpf_param_mem[3];
    LC3_INT  ltpf_mem_pitch;
    LC3_INT  ltpf_mem_pitch_fr;
    LC3_INT  ltpf_mem_beta_idx;
    LC3_INT  ltpf_mem_x_len;
    LC3_INT  ltpf_mem_y_len;
    LC3_INT  nbLostCmpt;
    LC3_INT  ltpf_conf_beta_idx;
    LC3_INT  spec_inv_idx;
    LC3_INT  concealMethod;
    LC3_INT  scf_idx[SCF_MAX_PARAM];
    uint8_t  resBits[MAX_RESBITS_LEN];
    LC3_INT  tns_idx[TNS_NUMFILTERS_MAX * MAXLAG];

    LC3_FLOAT ltpf_mem_x[3 * MAX_LEN];
    LC3_FLOAT ltpf_mem_y[3 * MAX_LEN];
    LC3_FLOAT ltpf_mem_gain;
    LC3_FLOAT plc_spec_prev[MAX_LEN];
    LC3_FLOAT ltpf_conf_beta;
    LC3_FLOAT plc_cum_alpha;
    LC3_FLOAT sqQdec_fl[MAX_LEN];
    LC3_FLOAT scf_q[M];
    LC3_FLOAT int_scf[MAX_BANDS_NUMBER];
    LC3_FLOAT x_fl[MAX_LEN];
    LC3_FLOAT imdct_mem[MAX_LEN];

    Dct4 dct4structImdct;
    LC3_INT plc_seed;
    
} DecSetup;

/* Constants and sampling rate derived values go in this struct */
struct LC3_Dec {
    DecSetup*  channel_setup[MAX_CHANNELS];
    const LC3_INT* W_fx;
    const LC3_INT* bands_offset;
    const LC3_INT* cutoffBins;

    LC3_INT fs;           /* sampling rate, 44.1 maps to 48 */
    LC3_INT fs_out;       /* output sampling rate */
    LC3_INT fs_idx;       /* sampling rate index */
    LC3_INT frame_length; /* sampling rate index */
    LC3_INT channels;     /* number of channels */
    LC3_FLOAT frame_ms;   /* frame length in ms (wrong for 44.1) */
    LC3_INT frame_dms;    /* frame length in ms * 10 (wrong for 44.1) */
    LC3_INT last_size;    /* size of last frame, without error protection */
    LC3_INT ep_enabled;   /* error protection enabled */
    LC3_INT error_report; /* corrected errors in last frame or -1 on error */

    LC3_INT          imdct_memLen;
    LC3_INT          imdct_winLen;
    LC3_INT          imdct_laZeros;
    const LC3_FLOAT* imdct_win;

    LC3_INT yLen;
    LC3_INT W_size;
    LC3_INT la_zeroes;
    LC3_INT bands_number;
    LC3_INT ltpf_mem_x_len;
    LC3_INT ltpf_mem_y_len;
    LC3_INT BW_cutoff_bits;
    LC3_INT plcMeth;
    LC3_INT tilt;
    
#ifdef ENABLE_HR_MODE
    LC3_INT hrmode;
    LC3_INT specflip;
#endif
    
};

#endif
