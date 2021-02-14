/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef SETUP_ENC_LC3_H
#define SETUP_ENC_LC3_H

#include "constants.h"

/* Channel state and bitrate-derived values go in this struct */
typedef struct
{
    Word16 *stEnc_mdct_mem; /* MDCT_MEM_LEN_MAX */
    Word32 *mdct_mem32;     /* MDCT_MEM_LEN_MAX */
    Word32  targetBitsOff;
    Word16  targetBytes;
    Word16  total_bits;
    Word16  targetBitsInit;
    Word16  targetBitsAri;
    Word16  enable_lpc_weighting;
    Word16  ltpf_enable;
    Word16  quantizedGainOff;
    Word16  tns_bits;
    Word16  targetBitsQuant;
    Word16  olpa_mem_s6k4_exp;
    Word16  olpa_mem_pitch;
    Word16  ltpf_mem_in_exp;
    Word16  ltpf_mem_normcorr;
    Word16  ltpf_mem_mem_normcorr;
    Word16  ltpf_mem_ltpf_on;
    Word16  ltpf_mem_pitch;
    Word16  mem_targetBits;
    Word16  mem_specBits;
    Word16  x_exp;
    Word16  resamp_exp;
    Word16  attack_handling; /* flag to enable attack handling */
    Word16  attdec_filter_mem[2];
    Word16  attdec_detected;
    Word16  attdec_position;
    Word32  attdec_acc_energy;
    Word16  attdec_scaling;
    Word32  resamp_mem32[60] ALIGN_BUFFER_STRUCT;
    Word32  r12k8_mem_50[2] ALIGN_BUFFER_STRUCT;
    Word16  r12k8_mem_in[60] ALIGN_BUFFER_STRUCT;
    Word16  r12k8_mem_out[24] ALIGN_BUFFER_STRUCT;
    Word16  olpa_mem_s12k8[3] ALIGN_BUFFER_STRUCT;
    Word16  olpa_mem_s6k4[LEN_6K4 + MAX_PITCH_6K4] ALIGN_BUFFER_STRUCT;
    Word16  ltpf_mem_in[LTPF_MEMIN_LEN + LEN_12K8 + 1] ALIGN_BUFFER_STRUCT;
    Word16 n_pccw;
    Word16 n_pc;
} EncSetup;

/* Constants and sampling rate derived values go in this struct */
struct LC3_Enc
{
    EncSetup *    channel_setup[MAX_CHANNELS];
    const Word16 *W_fx;
    const Word16 *bands_offset;

    Word32 fs;           /* encoder sampling rate 44.1 -> 48 */
    Word32 fs_in;        /* input sampling rate */
    Word32 bitrate;      /* global bitrate */
    Word16 fs_idx;       /* sampling rate index */
    Word16 frame_length; /* audio samples / frame */
    Word16 channels;     /* number of channels */
    Word16 epmode;       /* error protection mode */
    Word16 frame_dms;    /* frame length in dms (decimilliseconds, 10^-4)*/
    Word8  lc3_br_set;   /* indicate if bitrate has been set */

    Word16 yLen;
    Word16 W_size;
    Word16 la_zeroes;
    Word16 stEnc_mdct_mem_len;
    Word16 bands_number;
    Word16 nSubdivisions;
    Word16 ltpf_mem_in_len;
    Word16 envelope_bits;
    Word16 global_gain_bits;
    Word16 noise_fac_bits;
    Word16 BW_cutoff_bits;
    Word16 r12k8_mem_in_len;
    Word16 r12k8_mem_out_len;

    Word16 epmr;
    Word16 combined_channel_coding;
    Word32 bandwidth;
    Word16 bw_ctrl_cutoff_bin;
    Word16 bw_index;
};

#endif
