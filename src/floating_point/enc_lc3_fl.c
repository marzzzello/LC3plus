/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static void Enc_LC3_Channel_fl(LC3_Enc* encoder, int channel, int32_t* s_in, uint8_t* bytes, int bps)
{
    EncSetup* h_EncSetup;

    LC3_INT s_12k8_len = 0, T0_out = 0, ltpfBits = 0, BW_cutoff_idx = 0, tns_numfilters = 0, quantizedGain = 0,
        quantizedGainMin = 0, nbits = 0, nbits2 = 0, lastnz = 0, lsbMode = 0, gainChange = 0, bp_side = 0,
        mask_side = 0, fac_ns_idx = 0, numResBits = 0, tns_order[2] = {0}, i = 0;
    LC3_FLOAT normcorr = 0, gain = 0;

    LC3_FLOAT d_fl[MAX_LEN]                        = {0};
    LC3_INT   q_d[MAX_LEN]                         = {0};
    LC3_INT   indexes[TNS_NUMFILTERS_MAX * MAXLAG] = {0};

    h_EncSetup = encoder->channel_setup[channel];
    memset(bytes, 0, sizeof(uint8_t) * h_EncSetup->targetBytes);

    if (bps == 24) {
        for (i = 0; i < encoder->frame_length; i++) {
            h_EncSetup->s_in_scaled[i] = (LC3_FLOAT)(s_in[i] / LC3_POW(2, 8));
        }
    } else if (bps == 32) {
        for (i = 0; i < encoder->frame_length; i++) {
            h_EncSetup->s_in_scaled[i] = (LC3_FLOAT)(s_in[i] / LC3_POW(2, 16));
        }
    } else if (bps == 16) {
        for (i = 0; i < encoder->frame_length; i++) {
            h_EncSetup->s_in_scaled[i] = (LC3_FLOAT)((int16_t*)s_in)[i];
        }
    }

    /* MDCT */
    processMdct_fl(h_EncSetup->s_in_scaled, d_fl, &h_EncSetup->mdctStruct);
    
	/* 12.8 kHz resampler */
	process_resamp12k8_fl(h_EncSetup->s_in_scaled, encoder->frame_length, h_EncSetup->r12k8_mem_in,
						  encoder->r12k8_mem_in_len, h_EncSetup->r12k8_mem_50, h_EncSetup->r12k8_mem_out,
						  encoder->r12k8_mem_out_len, h_EncSetup->s_12k8, &s_12k8_len, encoder->fs_idx,
						  encoder->frame_dms, encoder->fs);

	/* Pitch estimation */
	processOlpa_fl(h_EncSetup->s_12k8, h_EncSetup->olpa_mem_s12k8, h_EncSetup->olpa_mem_s6k4,
				   &h_EncSetup->olpa_mem_pitch, &T0_out, &normcorr, s_12k8_len);

	/* LTPF encoder */
	process_ltpf_coder_fl(h_EncSetup->s_12k8, s_12k8_len + 1, h_EncSetup->ltpf_enable, T0_out, normcorr,
						  encoder->frame_dms, h_EncSetup->ltpf_mem_in, encoder->ltpf_mem_in_len,
						  &h_EncSetup->ltpf_mem_normcorr, &h_EncSetup->ltpf_mem_ltpf_on,
						  &h_EncSetup->ltpf_mem_pitch, h_EncSetup->ltpf_param, &h_EncSetup->ltpf_mem_mem_normcorr,
						  &ltpfBits);

    /* Attack detector */
    attack_detector_fl(h_EncSetup->s_in_scaled, encoder->frame_length, &h_EncSetup->attdec_position,
                       &h_EncSetup->attdec_acc_energy, &h_EncSetup->attdec_detected, h_EncSetup->attdec_filter_mem,
                       h_EncSetup->attack_handling);

    /* Per-band energy */
    processPerBandEnergy_fl(encoder->bands_number, encoder->bands_offset, h_EncSetup->ener, d_fl);
    
    /* Bandwidth cut-off detection */
#ifdef ENABLE_HR_MODE
    /* No BW Cutoff for 8 kHz and 96 kHz */
    if (encoder->fs_idx > 0 && encoder->hrmode == 0) {
#else
    if (encoder->fs_idx > 0) {
#endif
        processDetectCutoffWarped_fl(h_EncSetup->ener, encoder->fs_idx, encoder->frame_dms, &BW_cutoff_idx);
    } else {
        BW_cutoff_idx = encoder->fs_idx;
    }

    processSnsComputeScf_fl(h_EncSetup->ener, encoder->tilt, encoder->bands_number, h_EncSetup->scf,
                            h_EncSetup->attdec_detected, encoder->sns_damping);

    /* SNS Quantizer */
    process_snsQuantizesScf_Enc(h_EncSetup->scf, h_EncSetup->L_scf_idx, h_EncSetup->scf_q, h_EncSetup->dct2StructSNS);

    /* SNS Interpolation */
    processSnsInterpolateScf_fl(h_EncSetup->scf_q, 1, encoder->bands_number, h_EncSetup->int_scf);

    /* MDCT shaping */
    processMdctShaping_fl(d_fl, h_EncSetup->int_scf, encoder->bands_offset, encoder->bands_number);

    /* Bandwidth controller */
    if (encoder->bandwidth) {
        process_cutoff_bandwidth(d_fl, encoder->yLen, encoder->bw_ctrl_cutoff_bin);
        BW_cutoff_idx = MIN(BW_cutoff_idx, encoder->bw_index);
    }
        
    /* TNS encoder */
    processTnsCoder_fl(d_fl, BW_cutoff_idx, encoder->cutoffBins[BW_cutoff_idx], encoder->fs, encoder->frame_length,
                       encoder->frame_dms, h_EncSetup->total_bits, tns_order, indexes, &tns_numfilters,
                       &(h_EncSetup->tns_bits));

    /* Global Gain Estimation */
    h_EncSetup->targetBitsQuant = h_EncSetup->targetBitsInit - (h_EncSetup->tns_bits + ltpfBits);

    processEstimateGlobalGain_fl(d_fl, encoder->yLen, h_EncSetup->targetBitsQuant, &gain, &quantizedGain,
                                 &quantizedGainMin, h_EncSetup->quantizedGainOff, &h_EncSetup->targetBitsOff,
                                 &h_EncSetup->mem_targetBits, h_EncSetup->mem_specBits
#ifdef ENABLE_HR_MODE
                                 , encoder->hrmode, h_EncSetup->regBits, encoder->frame_ms
#endif

    );

    /* 1. Quantization */
    processQuantizeSpec_fl(d_fl, gain, q_d, encoder->yLen, h_EncSetup->total_bits, &nbits, &nbits2, encoder->fs,
                           &lastnz, h_EncSetup->codingdata, &lsbMode, -1, h_EncSetup->targetBitsQuant
#ifdef ENABLE_HR_MODE
		, encoder->hrmode
#endif
    );

    h_EncSetup->mem_specBits = nbits;

    /* Global Gain Adjustment */
    processAdjustGlobalGain_fl(&quantizedGain, quantizedGainMin, h_EncSetup->quantizedGainOff, &gain,
                               h_EncSetup->targetBitsQuant, h_EncSetup->mem_specBits, &gainChange, encoder->fs_idx);

    /* 2. Quantization */
    if (gainChange) {
        processQuantizeSpec_fl(d_fl, gain, q_d, encoder->yLen, h_EncSetup->total_bits, &nbits, &nbits2, encoder->fs,
                               &lastnz, h_EncSetup->codingdata, &lsbMode, 0, h_EncSetup->targetBitsQuant
#ifdef ENABLE_HR_MODE
		, encoder->hrmode
#endif
        );
    }

    /* Noise factor */
    processNoiseFactor_fl(&fac_ns_idx, d_fl, q_d, gain, encoder->cutoffBins[BW_cutoff_idx], encoder->frame_dms,
                          h_EncSetup->targetBytes
    );

    /* Residual Coding */
    if (lsbMode == 0) {
        processResidualCoding_fl(d_fl, q_d, gain, encoder->yLen, h_EncSetup->targetBitsQuant, nbits2,
                                 h_EncSetup->resBits, &numResBits
#ifdef ENABLE_HR_MODE
								 , encoder->hrmode
#endif
        );
    } else {
        numResBits = 0;
    }

    /* Entropy encoding */
    processEncoderEntropy_fl(bytes, &bp_side, &mask_side, h_EncSetup->targetBytes, encoder->BW_cutoff_bits,
                             BW_cutoff_idx, lastnz, encoder->yLen, lsbMode, quantizedGain, tns_numfilters, tns_order,
                             h_EncSetup->ltpf_param, h_EncSetup->L_scf_idx, fac_ns_idx
    );

    /* Artithmetic encoding */
    processAriEncoder_fl(bytes, bp_side, mask_side, q_d, tns_order, tns_numfilters, indexes, lastnz,
                         h_EncSetup->codingdata, h_EncSetup->resBits, numResBits, lsbMode, h_EncSetup->targetBitsAri,
                         h_EncSetup->enable_lpc_weighting);
}

int Enc_LC3_fl(LC3_Enc* encoder, void** input, uint8_t* output, int bps)
{
    int      ch = 0, output_size = 0;
    uint8_t* lc3buf = output;

    for (ch = 0; ch < encoder->channels; ch++) {
        Enc_LC3_Channel_fl(encoder, ch, input[ch], lc3buf, bps);
        lc3buf += encoder->channel_setup[ch]->targetBytes;
        output_size += encoder->channel_setup[ch]->targetBytes;
    }

    return output_size;
}
