/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static int Dec_LC3_Channel_fl(LC3_Dec* decoder, int channel, uint8_t* bs_in, void* s_out, int bps, int bfi_ext)
{
    DecSetup* h_DecSetup;
    LC3_INT       mask_side = 0, bp_side = 0, bfi = 0, gg_idx = 0, fac_ns_idx = 0, tns_numfilters = 0, bw_cutoff_idx = 0,
        lastnz = 0, lsbMode = 0, nf_seed = 0, zero_frame = 0, residualPresent = 0, nbits_residual = 0, bitsRead = 0,
        i = 0, tns_order[2] = {0}, sqQdec[MAX_LEN] = {0};

    h_DecSetup = decoder->channel_setup[channel];

    bfi = bfi_ext;

    /* Entropy decoding */
    if (!bfi) {
        processDecoderEntropy_fl(bs_in, h_DecSetup->targetBytes, &mask_side, &bp_side, decoder->yLen, decoder->fs_idx,
                                 decoder->BW_cutoff_bits, &bfi, &gg_idx, h_DecSetup->scf_idx, &fac_ns_idx,
                                 &tns_numfilters, tns_order, h_DecSetup->ltpf_param, &bw_cutoff_idx, &lastnz, &lsbMode, decoder->frame_dms
        );
    }
    
    /* Arithmetic decoding */
    if (!bfi) {
        processAriDecoder_fl(bs_in, bp_side, mask_side, decoder->yLen, decoder->fs_idx,
                             h_DecSetup->enable_lpc_weighting, tns_numfilters, lsbMode, lastnz, &bfi, tns_order,
                             fac_ns_idx, gg_idx, h_DecSetup->resBits, sqQdec, &nf_seed, h_DecSetup->tns_idx,
                             &zero_frame, h_DecSetup->targetBytes, &nbits_residual, &residualPresent
#ifdef ENABLE_HR_MODE
                             ,
                             decoder->hrmode
#endif
        );
        
        /* Cast from int to float */
        for (i = 0; i < decoder->yLen; i++) {
            h_DecSetup->sqQdec_fl[i] = (LC3_FLOAT)sqQdec[i];
        }
    }
    
    if (bfi != 1)
    {
        /* SNS Quantize Decoder */
        process_snsQuantizesScf_Dec(h_DecSetup->scf_idx, h_DecSetup->scf_q);
        
    }

    /* Decoding only if no bit error detected */
    if (!bfi) {
        /* Residual decoding */
        if (residualPresent) {
            processResidualDecoding_fl(&bitsRead, h_DecSetup->sqQdec_fl, decoder->yLen, h_DecSetup->resBits,
                                       nbits_residual
#ifdef ENABLE_HR_MODE
									   , decoder->hrmode
#endif
            );
        }
        
        /* Noise filling */
        if (zero_frame == 0) {
            processNoiseFilling_fl(h_DecSetup->sqQdec_fl, nf_seed, fac_ns_idx, decoder->cutoffBins[bw_cutoff_idx],
                                   decoder->frame_dms);
        }
        
        /* Application of global gain */
        processApplyGlobalGain_fl(h_DecSetup->sqQdec_fl, decoder->yLen, gg_idx, h_DecSetup->quantizedGainOff);

        /* TNS decoder */
        processTnsDecoder_fl(h_DecSetup->sqQdec_fl, h_DecSetup->tns_idx, tns_order, tns_numfilters,
                             decoder->cutoffBins[bw_cutoff_idx], decoder->frame_length, decoder->fs, decoder->frame_dms);

        /* SNS interpolation */
        processSnsInterpolateScf_fl(h_DecSetup->scf_q, 0, decoder->bands_number, h_DecSetup->int_scf);

        /* MDCT shaping */
        processMdctShaping_fl(h_DecSetup->sqQdec_fl, h_DecSetup->int_scf, decoder->bands_offset, decoder->bands_number);
    }
    
    
    
    if (bfi) {
        /* Apply PLC */
        processPlcApply(h_DecSetup->plc_spec_prev, decoder->yLen, &h_DecSetup->nbLostCmpt, &h_DecSetup->plc_cum_alpha,
                        &h_DecSetup->plc_seed, h_DecSetup->sqQdec_fl);
    } else {
        /* Update PLC */
        processPlcUpdate(h_DecSetup->sqQdec_fl, h_DecSetup->plc_spec_prev, decoder->yLen, &h_DecSetup->nbLostCmpt,
                         &h_DecSetup->plc_cum_alpha);
    }
    
    /* IMDCT */
    {
        ProcessingIMDCT_fl(h_DecSetup->sqQdec_fl, decoder->frame_length, decoder->imdct_win, decoder->imdct_winLen, decoder->imdct_laZeros, h_DecSetup->imdct_mem, h_DecSetup->x_fl, &h_DecSetup->dct4structImdct);
    }
    
    
    /* LTPF decoder */
    if (decoder->fs_idx != 5)
    {
        process_ltpf_decoder_fl(h_DecSetup->x_fl, decoder->frame_length, h_DecSetup->x_fl, decoder->fs,
                                h_DecSetup->ltpf_mem_x, h_DecSetup->ltpf_mem_y, &h_DecSetup->ltpf_mem_pitch,
                                &h_DecSetup->ltpf_mem_pitch_fr, &h_DecSetup->ltpf_mem_gain, &h_DecSetup->ltpf_mem_beta_idx,
                                bfi, h_DecSetup->ltpf_param, h_DecSetup->ltpf_param_mem, h_DecSetup->ltpf_conf_beta_idx,
                                h_DecSetup->ltpf_conf_beta,
                                h_DecSetup->concealMethod, 
                                h_DecSetup->plc_cum_alpha);
    } else {
        memmove(h_DecSetup->x_fl, h_DecSetup->x_fl, decoder->frame_length * sizeof(LC3_FLOAT));
    }

    {
        /* Round, scale and copy output to output buffer */
        if (bps == 16) {
            for (i = 0; i < decoder->frame_length; i++) {
                LC3_FLOAT tmp            = round(LC3_POW(2, 16 - 1) * (h_DecSetup->x_fl[i] * LC3_POW(2, -15)));
                ((int16_t*)s_out)[i] = (int16_t)fmaxf(fminf(tmp, 32767), -32768);
            }
        } else {
            for (i = 0; i < decoder->frame_length; i++) {
                ((int32_t*)s_out)[i] = (int32_t)round(LC3_POW(2, bps - 1) * (h_DecSetup->x_fl[i] * LC3_POW(2, -15)));
            }
        }
    }

    return bfi;
}

/* num_bytes = 0 -> bad frame */
LC3_Error Dec_LC3_fl(LC3_Dec* decoder, uint8_t* input, int num_bytes, void** output, int bps)
{
    int       ch = 0, bfi = !num_bytes;
    LC3_Error err = LC3_OK;

    for (ch = 0; ch < decoder->channels; ch++) {
#ifdef ENABLE_PADDING
        if (!bfi) {
            LC3_INT padding_len;

            bfi = paddingDec_fl(input, num_bytes, decoder->yLen, decoder->BW_cutoff_bits, &padding_len);

            num_bytes = num_bytes - padding_len;
            if (num_bytes < 20 || num_bytes > MAX_NBYTES) {
                bfi = 1;  /* frame below the minimum, frame is broken */
            }
        }
#endif

        if (!bfi && num_bytes != decoder->last_size) {
            err = update_dec_bitrate(decoder, ch, num_bytes);
            if (err)
                return err;
            decoder->last_size = num_bytes;
        }

        bfi = Dec_LC3_Channel_fl(decoder, ch, input, output[ch], bps, bfi);
        input += decoder->channel_setup[ch]->targetBytes;
    }

    return bfi ? LC3_DECODE_ERROR : LC3_OK;
}

