/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef SETUP_ENC_LC3_FL_H
#define SETUP_ENC_LC3_FL_H

#include "constants.h"

/* Channel state and bitrate-derived values go in this struct */
typedef struct {
    LC3_FLOAT targetBitsOff;
    LC3_FLOAT ltpf_mem_normcorr;
    LC3_FLOAT ltpf_mem_mem_normcorr;
    LC3_FLOAT attdec_filter_mem[2];
    LC3_FLOAT attdec_acc_energy;
    LC3_FLOAT r12k8_mem_50[2];
    LC3_FLOAT r12k8_mem_in[120];
    LC3_FLOAT r12k8_mem_out[24];
    LC3_FLOAT olpa_mem_s12k8[3];
    LC3_FLOAT olpa_mem_s6k4[LEN_6K4 + MAX_PITCH_6K4];
    LC3_FLOAT ltpf_mem_in[LTPF_MEMIN_LEN + LEN_12K8 + 1];
    LC3_FLOAT s_in_scaled[MAX_LEN];
    LC3_FLOAT s_12k8[LEN_12K8 + 1];
    LC3_FLOAT ener[MAX_BANDS_NUMBER];
    LC3_FLOAT scf_q[M];
    LC3_FLOAT scf[M];
    LC3_FLOAT int_scf[MAX_BANDS_NUMBER];

    LC3_INT targetBytes;
    LC3_INT total_bits;
    LC3_INT targetBitsInit;
    LC3_INT targetBitsAri;
    LC3_INT enable_lpc_weighting;
    LC3_INT ltpf_enable;
    LC3_INT quantizedGainOff;
    LC3_INT tns_bits;
    LC3_INT targetBitsQuant;
    LC3_INT olpa_mem_pitch;
    LC3_INT ltpf_mem_ltpf_on;
    LC3_INT ltpf_mem_pitch;
    LC3_INT mem_targetBits;
    LC3_INT mem_specBits;
    LC3_INT attack_handling; /* flag to enable attack handling */
    LC3_INT attdec_detected;
    LC3_INT attdec_position;
    LC3_INT ltpf_param[3];
    LC3_INT L_scf_idx[SCF_MAX_PARAM];
    LC3_INT codingdata[3 * MAX_LEN];
#ifdef ENABLE_HR_MODE
    uint8_t resBits[MAX_RESBITS_LEN];
    LC3_INT regBits;
#else
    uint8_t resBits[MAX_RESBITS_LEN];
#endif

    Mdct mdctStruct;
    Dct2 dct2StructSNS;
} EncSetup;

/* Constants and sampling rate derived values go in this struct */
struct LC3_Enc {
    EncSetup*  channel_setup[MAX_CHANNELS];
    const LC3_INT* W_fx;
    const LC3_INT* bands_offset;
    const LC3_INT* cutoffBins;

    LC3_INT fs;           /* encoder sampling rate 44.1 -> 48 */
    LC3_INT fs_in;        /* input sampling rate */
    LC3_INT bitrate;      /* global bitrate */
    LC3_INT fs_idx;       /* sampling rate index */
    LC3_INT frame_length; /* audio samples / frame */
    LC3_INT channels;     /* number of channels */
    LC3_INT epmode;      /* error protection mode */
    LC3_FLOAT frame_ms;   /* frame length in ms (wrong for 44.1) */
    LC3_INT frame_dms;    /* frame length in ms * 10 (wrong for 44.1) */
    LC3_INT tilt;
    LC3_INT lc3_br_set;
    LC3_INT yLen;
    LC3_INT W_size;
    LC3_INT la_zeroes;
    LC3_INT stEnc_mdct_mem_len;
    LC3_INT bands_number;
    LC3_INT nSubdivisions;
    LC3_INT ltpf_mem_in_len;
    LC3_INT envelope_bits;
    LC3_INT global_gain_bits;
    LC3_INT noise_fac_bits;
    LC3_INT BW_cutoff_bits;
    LC3_INT r12k8_mem_in_len;
    LC3_INT r12k8_mem_out_len;
#ifdef ENABLE_HR_MODE
    LC3_INT hrmode;
#endif
    LC3_INT bandwidth;
    LC3_INT bw_ctrl_cutoff_bin;
    LC3_INT bw_index;
    LC3_FLOAT sns_damping;
};

#endif
