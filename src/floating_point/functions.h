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

#include "clib.h"
#include "defines.h"
#include "float.h"
#include "lc3.h"
#include "setup_dec_lc3.h"
#include "setup_enc_lc3.h"
#include "structs.h"
#include "util.h"


/* FFT */
#include "fft/iis_fft.h"
#include "fft/iisfft.h"

/* fft.c */
void fft_init(Fft* fft, LC3_INT length);
void fft_free(Fft* fft);
void fft_apply(Fft* fft, const Complex* input, Complex* output);

/* dct.c */
void dct2_init(Dct2* dct, LC3_INT length);
void dct2_free(Dct2* dct);
void dct2_apply(Dct2* dct, const LC3_FLOAT* input, LC3_FLOAT* output);

void dct3_init(Dct3* dct, LC3_INT length);
void dct3_free(Dct3* dct);
void dct3_apply(Dct3* dct, const LC3_FLOAT* input, LC3_FLOAT* output);

void dct4_init(Dct4* dct, LC3_INT length);
void dct4_free(Dct4* dct);
void dct4_apply(Dct4* dct, const LC3_FLOAT* input, LC3_FLOAT* output);

/* mdct.c */
#ifdef ENABLE_HR_MODE
void mdct_init(Mdct* mdct, LC3_INT length, LC3_INT frame_dms, LC3_INT hrmode);
#else
void mdct_init(Mdct* mdct, LC3_INT length, LC3_INT frame_dms);
#endif
void mdct_free(Mdct* mdct);
void mdct_apply(const LC3_FLOAT* input, LC3_FLOAT* output, Mdct* mdct);

#ifdef ENABLE_PADDING
LC3_INT paddingDec_fl(LC3_UINT8* bytes, LC3_INT numbytes, LC3_INT N, LC3_INT bw_cutoff_bits, LC3_INT* total_padding);
#endif

#ifdef ENABLE_HR_MODE

void processEncoderEntropy_fl(LC3_UINT8* bytes, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT numbytes, LC3_INT bw_cutoff_bits,
                              LC3_INT bw_cutoff_idx, LC3_INT lastnz, LC3_INT N, LC3_INT lsbMode, LC3_INT gg_idx, LC3_INT num_tns_filters,
                              LC3_INT* tns_order, LC3_INT* ltpf_idx, LC3_INT* scf_idx, LC3_INT fac_ns_idx
							  );
void processDecoderEntropy_fl(LC3_UINT8* bytes, LC3_INT numbytes, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT N, LC3_INT fs_idx,
                              LC3_INT bw_cutoff_bits, LC3_INT* bfi, LC3_INT* gg_idx, LC3_INT* scf_idx, LC3_INT* fac_ns_idx,
                              LC3_INT* tns_numfilters, LC3_INT* tns_order, LC3_INT* ltpf_idx, LC3_INT* bw_cutoff_idx, LC3_INT* lastnz,
                              LC3_INT* lsbMode, LC3_INT frame_dms
							  );
void processQuantizeSpec_fl(LC3_FLOAT x[], LC3_FLOAT gain, LC3_INT xq[], LC3_INT nt, LC3_INT totalBits, LC3_INT* nbits, LC3_INT* nbits2, LC3_INT fs,
                            LC3_INT* lastnzout, LC3_INT* codingdata, LC3_INT* lsbMode, LC3_INT mode, LC3_INT target, LC3_INT hrmode
							);
void processEstimateGlobalGain_fl(LC3_FLOAT x[], LC3_INT lg, LC3_INT nbitsSQ, LC3_FLOAT* gain, LC3_INT* quantizedGain,
                                  LC3_INT* quantizedGainMin, LC3_INT quantizedGainOff, LC3_FLOAT* targetBitsOff,
                                  LC3_INT* old_targetBits, LC3_INT old_specBits, LC3_INT bq_mode
#ifdef ENABLE_HR_MODE
								  , LC3_INT regBits, LC3_FLOAT frame_ms
#endif
);

void processAriDecoder_fl(LC3_UINT8* bytes, LC3_INT bp_side, LC3_INT mask_side, LC3_INT L_spec, LC3_INT fs_idx, LC3_INT enable_lpc_weighting,
                          LC3_INT tns_numfilters, LC3_INT lsbMode, LC3_INT lastnz, LC3_INT* bfi, LC3_INT* tns_order, LC3_INT fac_ns_idx,
                          LC3_INT gg_idx, uint8_t* resBits, LC3_INT* x, LC3_INT* nf_seed, LC3_INT* tns_idx, LC3_INT* zero_frame, LC3_INT numbytes,
                          LC3_INT* nbits_residual, LC3_INT* residualPresent, LC3_INT bq_mode);
#else

void processEncoderEntropy_fl(LC3_UINT8* bytes, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT numbytes, LC3_INT bw_cutoff_bits,
                              LC3_INT bw_cutoff_idx, LC3_INT lastnz, LC3_INT N, LC3_INT lsbMode, LC3_INT gg_idx, LC3_INT num_tns_filters,
                              LC3_INT* tns_order, LC3_INT* ltpf_idx, LC3_INT* scf_idx, LC3_INT fac_ns_idx);
void processDecoderEntropy_fl(LC3_UINT8* bytes, LC3_INT numbytes, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT N, LC3_INT fs_idx,
                              LC3_INT bw_cutoff_bits, LC3_INT* bfi, LC3_INT* gg_idx, LC3_INT* scf_idx, LC3_INT* fac_ns_idx,
                              LC3_INT* tns_numfilters, LC3_INT* tns_order, LC3_INT* ltpf_idx, LC3_INT* bw_cutoff_idx, LC3_INT* lastnz,
                              LC3_INT* lsbMode, LC3_INT frame_dms);
void processQuantizeSpec_fl(LC3_FLOAT x[], LC3_FLOAT gain, LC3_INT xq[], LC3_INT nt, LC3_INT totalBits, LC3_INT* nbits, LC3_INT* nbits2, LC3_INT fs,
                            LC3_INT* lastnzout, LC3_INT* codingdata, LC3_INT* lsbMode, LC3_INT mode, LC3_INT target);
void processEstimateGlobalGain_fl(LC3_FLOAT x[], LC3_INT lg, LC3_INT nbitsSQ, LC3_FLOAT* gain, LC3_INT* quantizedGain,
                                  LC3_INT* quantizedGainMin, LC3_INT quantizedGainOff, LC3_FLOAT* targetBitsOff,
                                  LC3_INT* old_targetBits, LC3_INT old_specBits);
void processAriDecoder_fl(LC3_UINT8* bytes, LC3_INT bp_side, LC3_INT mask_side, LC3_INT L_spec, LC3_INT fs_idx, LC3_INT enable_lpc_weighting,
                          LC3_INT tns_numfilters, LC3_INT lsbMode, LC3_INT lastnz, LC3_INT* bfi, LC3_INT* tns_order, LC3_INT fac_ns_idx,
                          LC3_INT gg_idx, uint8_t* resBits, LC3_INT* x, LC3_INT* nf_seed, LC3_INT* tns_idx, LC3_INT* zero_frame, LC3_INT numbytes,
                          LC3_INT* nbits_residual, LC3_INT* residualPresent);
#endif

void processMdctShaping_fl(LC3_FLOAT x[], LC3_FLOAT gains[], const LC3_INT bands_offset[], LC3_INT fdns_npts);

void processResidualCoding_fl(LC3_FLOAT x[], LC3_INT xq[], LC3_FLOAT gain, LC3_INT L_spec, LC3_INT targetBits, LC3_INT nBits, uint8_t * resBits,
                              LC3_INT* numResBits
#ifdef ENABLE_HR_MODE
							  , LC3_INT hrmode
#endif
);

void processResidualDecoding_fl(LC3_INT* bitsRead, LC3_FLOAT x[], LC3_INT L_spec, uint8_t prm[], LC3_INT resQBits
#ifdef ENABLE_HR_MODE
								, LC3_INT hrmode
#endif
);

void processAdjustGlobalGain_fl(LC3_INT* gg_idx, LC3_INT gg_idx_min, LC3_INT gg_idx_off, LC3_FLOAT* gain, LC3_INT target, LC3_INT nBits,
                                LC3_INT* gainChange, LC3_INT fs_idx);

void processApplyGlobalGain_fl(LC3_FLOAT x[], LC3_INT xLen, LC3_INT global_gain_idx, LC3_INT global_gain_off);

void processNoiseFactor_fl(LC3_INT* fac_ns_idx, LC3_FLOAT x[], LC3_INT xq[], LC3_FLOAT gg, LC3_INT BW_cutoff_idx, LC3_INT frame_dms,
                           LC3_INT target_bytes
							);

void processNoiseFilling_fl(LC3_FLOAT xq[], LC3_INT nfseed, LC3_INT fac_ns_idx, LC3_INT BW_cutoff_idx, LC3_INT frame_dmss);

void processOlpa_fl(LC3_FLOAT* wsp, LC3_FLOAT* mem_lp_decim2, LC3_FLOAT* mem_old_d_wsp, LC3_INT* mem_old_T0, LC3_INT* T0_out,
                    LC3_FLOAT* normcorr_out, LC3_INT wsp_len);

void processTnsCoder_fl(LC3_FLOAT* x, LC3_INT bw_cutoff_idx, LC3_INT bw_fcbin, LC3_INT fs, LC3_INT N, LC3_INT frame_dms, LC3_INT nBits,
                        LC3_INT* order_out, LC3_INT* rc_idx, LC3_INT* tns_numfilters, LC3_INT* bits_out);
void levinsonDurbin(LC3_FLOAT* r, LC3_FLOAT* out_lev, LC3_FLOAT* rc_unq, LC3_FLOAT* error, LC3_INT len);

void processTnsDecoder_fl(LC3_FLOAT* x, LC3_INT* rc_idx, LC3_INT* order, LC3_INT numfilters, LC3_INT bw_fcbin, LC3_INT N, LC3_INT fs, LC3_INT frame_dms);

void processSnsComputeScf_fl(LC3_FLOAT* x, LC3_INT tilt, LC3_INT xLen, LC3_FLOAT* gains, LC3_INT smooth, LC3_FLOAT sns_damping);

void processSnsInterpolateScf_fl(LC3_FLOAT* gains, LC3_INT encoder_side, LC3_INT bands_number, LC3_FLOAT* gains_LC3_INT);

void processDetectCutoffWarped_fl(LC3_FLOAT* d2, LC3_INT fs_idx, LC3_INT frame_dms, LC3_INT* bw_idx);

void processPerBandEnergy_fl(LC3_INT bands_number, const LC3_INT* acc_coeff_per_band, LC3_FLOAT* d2, LC3_FLOAT* d);

void ProcessingIMDCT_fl(LC3_FLOAT* y, LC3_INT yLen, const LC3_FLOAT* win, LC3_INT winLen, LC3_INT last_zeros, LC3_FLOAT* mem, LC3_FLOAT* x,
                        Dct4* dct);

void process_ltpf_coder_fl(LC3_FLOAT* x, LC3_INT xLen, LC3_INT ltpf_enable, LC3_INT pitch_ol, LC3_FLOAT pitch_ol_norm_corr, LC3_INT frame_dms,
                           LC3_FLOAT* mem_old_x, LC3_INT memLen, LC3_FLOAT* mem_norm_corr_past, LC3_INT* mem_on, LC3_INT* mem_pitch,
                           LC3_INT* param, LC3_FLOAT* mem_norm_corr_past_past, LC3_INT* bits);

void process_ltpf_decoder_fl(LC3_FLOAT* x, LC3_INT xLen, LC3_FLOAT* y, LC3_INT fs, LC3_FLOAT* mem_old_x, LC3_FLOAT* mem_old_y,
                             LC3_INT* mem_pitch_LC3_INT, LC3_INT* mem_pitch_fr, LC3_FLOAT* mem_gain, LC3_INT* mem_beta_idx, LC3_INT bfi,
                             LC3_INT* param, LC3_INT* mem_param, LC3_INT conf_beta_idx, LC3_FLOAT conf_beta, LC3_INT concealMethod, LC3_FLOAT damping);

void process_resamp12k8_fl(LC3_FLOAT x[], LC3_INT x_len, LC3_FLOAT mem_in[], LC3_INT mem_in_len, LC3_FLOAT mem_50[], LC3_FLOAT mem_out[],
                           LC3_INT mem_out_len, LC3_FLOAT y[], LC3_INT* y_len, LC3_INT fs_idx, LC3_INT frame_dms, LC3_INT fs);

void write_bit_backward_fl(LC3_UINT8* ptr, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT bit);
void write_uint_backward_fl(LC3_UINT8* ptr, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT val, LC3_INT numbits);


void processAriEncoder_fl(LC3_UINT8* bytes, LC3_INT bp_side, LC3_INT mask_side, LC3_INT* x, LC3_INT* tns_order, LC3_INT tns_numfilters,
                          LC3_INT* tns_idx, LC3_INT lastnz, LC3_INT* codingdata, uint8_t* res_bits, LC3_INT resBitsLen, LC3_INT lsbMode,
                          LC3_INT nbbits, LC3_INT enable_lpc_weighting);

void attack_detector_fl(LC3_FLOAT* in, LC3_INT frame_size, LC3_INT* lastAttackPosition, LC3_FLOAT* accNrg, LC3_INT* attackFlag,
                        LC3_FLOAT* attdec_filter_mem, LC3_INT attackHandlingOn);

void process_snsQuantizesScf_Enc(LC3_FLOAT* env, LC3_INT* index, LC3_FLOAT* envq, Dct2 dct2structSNS);

void process_snsQuantizesScf_Dec(LC3_INT* scf_idx, LC3_FLOAT* scf_q);

void processMdct_fl(LC3_FLOAT* in, LC3_FLOAT* out, Mdct* mdctStruct);

int       alloc_encoder(LC3_Enc* encoder, int channels);
void      set_enc_frame_params(LC3_Enc* encoder);
LC3_Error update_enc_bitrate(LC3_Enc* encoder, int bitrate);

LC3_Error FillEncSetup(LC3_Enc* encoder, int samplerate, int channels);

/* Setup Functions */
int       alloc_decoder(LC3_Dec* decoder, int channels);
void      set_dec_frame_params(LC3_Dec* decoder);
LC3_Error update_dec_bitrate(LC3_Dec* decoder, int ch, int nBytes);
LC3_Error FillDecSetup(LC3_Dec* decoder, int samplerate, int channels, LC3_PlcMode plc_mode);

int       Enc_LC3_fl(LC3_Enc* encoder, void** input, LC3_UINT8* output, int bps);
LC3_Error Dec_LC3_fl(LC3_Dec* decoder, LC3_UINT8* input, int input_bytes, void** output, int bps);

void* balloc(void* base, size_t* base_size, size_t size);




void processPlcApply(LC3_FLOAT* spec_prev, LC3_INT L_spec, LC3_INT* nbLostCmpt, LC3_FLOAT* cum_alpha, LC3_INT* seed, LC3_FLOAT* spec_out);
void processPlcUpdate(LC3_FLOAT* spec_cur, LC3_FLOAT* spec_prev, LC3_INT len, LC3_INT* nbLostCmpt, LC3_FLOAT* cum_alpha);

void process_cutoff_bandwidth(LC3_FLOAT* d_fl, LC3_INT len, LC3_INT bw_bin);
void update_enc_bandwidth(LC3_Enc* encoder, LC3_INT bandwidth);



#endif
