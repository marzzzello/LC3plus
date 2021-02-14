/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "setup_enc_lc3.h"
#include "functions.h"
#include <stdio.h>

/* if encoder is null only size is reported */
int alloc_encoder(LC3_Enc* encoder, int channels)
{
    int    ch   = 0;
    size_t size = sizeof(LC3_Enc);

    for (ch = 0; ch < channels; ch++) {
        EncSetup* setup = balloc(encoder, &size, sizeof(EncSetup));
        if (encoder) {
            encoder->channel_setup[ch] = setup;
        }
    }

    return (int)size;
}

LC3_Error FillEncSetup(LC3_Enc* encoder, int samplerate, int channels)
{
    memset(encoder, 0, lc3_enc_get_size(samplerate, channels));
    alloc_encoder(encoder, channels);

    encoder->fs     = CODEC_FS(samplerate);
    encoder->fs_in  = samplerate;
    encoder->fs_idx = FS2FS_IDX(encoder->fs);
    encoder->frame_dms = 100;
#ifdef ENABLE_HR_MODE
    if (encoder->fs_idx > 4) {
        encoder->fs_idx = 5;
    }
#endif
    encoder->channels          = channels;
    encoder->frame_ms          = 10;
    encoder->envelope_bits     = 38;
    encoder->global_gain_bits  = 8;
    encoder->noise_fac_bits    = 3;
    encoder->BW_cutoff_bits    = BW_cutoff_bits_all[encoder->fs_idx];

    encoder->r12k8_mem_in_len  = 2 * 8 * encoder->fs / 12800;
    encoder->r12k8_mem_out_len = 24;

    if (encoder->fs == 8000) {
        encoder->tilt = 14;
    } else if (encoder->fs == 16000) {
        encoder->tilt = 18;
    } else if (encoder->fs == 24000) {
        encoder->tilt = 22;
    } else if (encoder->fs == 32000) {
        encoder->tilt = 26;
    } else if (encoder->fs == 48000) {
        encoder->tilt = 30;
    }
#ifdef ENABLE_HR_MODE
    else if (encoder->fs == 96000) {
        encoder->tilt = 34;
    }
#endif

    set_enc_frame_params(encoder);
    return LC3_OK;
}

/* set frame config params */
void set_enc_frame_params(LC3_Enc* encoder)
{
    int       ch = 0;
    EncSetup* setup;

    encoder->frame_length       = ceil(encoder->fs * 10 / 1000); /* fs * 0.01*2^6 */
#ifdef ENABLE_HR_MODE
    if (encoder->hrmode == 1)
    {
    	encoder->yLen = encoder->frame_length;
        encoder->sns_damping = 0.6;
    }
    else
    {
    	encoder->yLen = MIN(MAX_BW, encoder->frame_length);
        encoder->sns_damping = 0.85;
    }
#else
    encoder->yLen               = MIN(MAX_BW, encoder->frame_length);
    encoder->sns_damping        = 0.85;
#endif
    encoder->stEnc_mdct_mem_len = encoder->frame_length - encoder->la_zeroes;
    encoder->bands_number       = 64;
    encoder->nSubdivisions      = 3;
    encoder->ltpf_mem_in_len    = LTPF_MEMIN_LEN;

#ifdef ENABLE_HR_MODE
    if (encoder->hrmode)
    {
    	encoder->BW_cutoff_bits    = 0;
    }
    else
    {
    	encoder->BW_cutoff_bits    = BW_cutoff_bits_all[encoder->fs_idx];
    }
#endif

    if (encoder->frame_ms == 10) {
#ifdef ENABLE_HR_MODE
        if (encoder->hrmode)
        {
            encoder->bands_offset = ACC_COEFF_PER_BAND_HR[encoder->fs_idx];
            encoder->cutoffBins   = BW_cutoff_bin_all_HR;
        } else
        {
            encoder->bands_offset = ACC_COEFF_PER_BAND[encoder->fs_idx];
            encoder->cutoffBins   = BW_cutoff_bin_all;
        }
#else
        encoder->bands_offset = ACC_COEFF_PER_BAND[encoder->fs_idx];
        encoder->cutoffBins   = BW_cutoff_bin_all;
#endif
    }
    else if (encoder->frame_ms == 2.5) {
#ifdef ENABLE_HR_MODE
        if (encoder->hrmode)
        {
            encoder->bands_offset = ACC_COEFF_PER_BAND_2_5ms_HR[encoder->fs_idx];
            encoder->cutoffBins   = BW_cutoff_bin_all_2_5ms_HR; 
        } else
        {
            encoder->bands_offset = ACC_COEFF_PER_BAND_2_5ms[encoder->fs_idx];
            encoder->cutoffBins   = BW_cutoff_bin_all_2_5ms;
        }
#else
        encoder->bands_offset = ACC_COEFF_PER_BAND_2_5ms[encoder->fs_idx];
        encoder->cutoffBins   = BW_cutoff_bin_all_2_5ms;
#endif
    }
    else if (encoder->frame_ms == 5) {
#ifdef ENABLE_HR_MODE
        if (encoder->hrmode)
        {
            encoder->bands_offset = ACC_COEFF_PER_BAND_5ms_HR[encoder->fs_idx];
            encoder->cutoffBins   = BW_cutoff_bin_all_5ms_HR;
        } else {
            encoder->bands_offset = ACC_COEFF_PER_BAND_5ms[encoder->fs_idx];
            encoder->cutoffBins   = BW_cutoff_bin_all_5ms;
        }
#else
        encoder->bands_offset = ACC_COEFF_PER_BAND_5ms[encoder->fs_idx];
        encoder->cutoffBins   = BW_cutoff_bin_all_5ms;
#endif
    }

    if (encoder->frame_ms == 2.5) {
        encoder->frame_length = encoder->frame_length >> 2;
#ifdef  ENABLE_HR_MODE
        encoder->yLen /= 4;
#else
        encoder->yLen = MIN(MAX_BW >> 2, encoder->frame_length);
#endif
        encoder->stEnc_mdct_mem_len = encoder->frame_length - encoder->la_zeroes;
#ifdef  ENABLE_HR_MODE
        if (encoder->hrmode)
        {
            encoder->bands_number       = bands_number_2_5ms_HR[encoder->fs_idx];
        } else
#endif
        {
            encoder->bands_number       = bands_number_2_5ms[encoder->fs_idx];
        }

        encoder->nSubdivisions      = 2;
        encoder->ltpf_mem_in_len    = LTPF_MEMIN_LEN + (LEN_12K8 >> 2) * 3;
    }
    
    
    if (encoder->frame_ms == 5) {
        encoder->frame_length = encoder->frame_length >> 1;
#ifdef  ENABLE_HR_MODE
        encoder->yLen /= 2;
#else
        encoder->yLen = MIN(MAX_BW >> 1, encoder->frame_length);
#endif
        encoder->stEnc_mdct_mem_len = encoder->frame_length - encoder->la_zeroes;
        encoder->bands_number       = bands_number_5ms[encoder->fs_idx];
        encoder->nSubdivisions      = 2;
        encoder->ltpf_mem_in_len    = LTPF_MEMIN_LEN + (LEN_12K8 >> 1);
    }

    for (ch = 0; ch < encoder->channels; ch++) {
        setup = encoder->channel_setup[ch];

        if (setup->mdctStruct.mem != NULL) {
            mdct_free(&setup->mdctStruct);
#ifdef ENABLE_HR_MODE
            mdct_init(&setup->mdctStruct, encoder->frame_length, encoder->frame_dms, encoder->hrmode);
#else
            mdct_init(&setup->mdctStruct, encoder->frame_length, encoder->frame_dms);
#endif

            dct2_free(&setup->dct2StructSNS);
            dct2_init(&setup->dct2StructSNS, M);
        } else {
#ifdef ENABLE_HR_MODE
        	mdct_init(&setup->mdctStruct, encoder->frame_length, encoder->frame_dms, encoder->hrmode);
#else
            mdct_init(&setup->mdctStruct, encoder->frame_length, encoder->frame_dms);
#endif
            dct2_init(&setup->dct2StructSNS, M);
        }
    }
}

/* change encoder bitrate */
LC3_Error update_enc_bitrate(LC3_Enc* encoder, int bitrate)
{
    int ch = 0, bitsTmp = 0, minBR = 0, maxBR = 0, nbytes = 0, channel_bytes;

#ifdef ENABLE_HR_MODE
    if (encoder->hrmode)
    {
        switch (encoder->frame_dms)
        {
        case 25:
            maxBR = 672000;
            if (encoder->fs == 48000) {minBR = 172800;}
            else if (encoder->fs == 96000) {minBR = 198400;}
            else { return LC3_HRMODE_ERROR;}
            break;
        case 50:
            maxBR = 600000;
            if (encoder->fs == 48000) {minBR = 148800;}
            else if (encoder->fs == 96000) {minBR = 174400;}
            else { return LC3_HRMODE_ERROR;}
            break;
        case 100:
            maxBR = 500000;
            if (encoder->fs == 48000) {minBR = 124800;}
            else if (encoder->fs == 96000) {minBR = 149600;}
            else { return LC3_HRMODE_ERROR;}
            break;
        default:
            return LC3_HRMODE_ERROR;
        }
    }
    else
#endif
    {
        minBR = MIN_NBYTES * 8 * (1000 / encoder->frame_ms);
        maxBR = (encoder->fs_in == 44100 ? MAX_NBYTES : MAX_NBYTES_RED) * 8 * (1000 / encoder->frame_ms);
    }
    minBR *= encoder->channels;
    maxBR *= encoder->channels;

    if (bitrate < minBR || bitrate > maxBR) {
        return LC3_BITRATE_ERROR;
    }
    
    encoder->lc3_br_set = 1;

    nbytes = bitrate * encoder->frame_length / (8 * encoder->fs_in);

    for (ch = 0; ch < encoder->channels; ch++) {

        EncSetup* setup = encoder->channel_setup[ch];
        channel_bytes = nbytes / encoder->channels;
        
        if(ch < nbytes % encoder->channels)
    	    {
            channel_bytes ++;
    	    }

        setup->targetBytes    = channel_bytes;
        setup->total_bits     = setup->targetBytes << 3;
        setup->targetBitsInit = setup->total_bits - encoder->envelope_bits - encoder->global_gain_bits -
                                encoder->noise_fac_bits - encoder->BW_cutoff_bits -
                                ceil(LC3_LOG2(encoder->frame_length / 2)) - 2 - 1;

        if (setup->total_bits > 1280) {
            setup->targetBitsInit = setup->targetBitsInit - 1;
        }
        if (setup->total_bits > 2560) {
            setup->targetBitsInit = setup->targetBitsInit - 1;
        }

#ifdef ENABLE_HR_MODE
        if (encoder->hrmode)
        {
        		setup->targetBitsInit -= 1;
        }
#endif

        setup->targetBitsAri        = setup->total_bits;
        setup->enable_lpc_weighting = setup->total_bits < 480;

        if (encoder->frame_ms == 5) {
            setup->enable_lpc_weighting = setup->total_bits < 240;
        }
        if (encoder->frame_ms == 2.5) {
            setup->enable_lpc_weighting = setup->total_bits < 120;
        }

        setup->quantizedGainOff =
            -(MIN(115, setup->total_bits / (10 * (encoder->fs_idx + 1))) + 105 + 5 * (encoder->fs_idx + 1));
        if (encoder->frame_ms == 10 && ((encoder->fs_in >= 44100 && setup->targetBytes >= 100) ||
                                        (encoder->fs_in == 32000 && setup->targetBytes >= 81)) && setup->targetBytes < 340) {
            setup->attack_handling = 1;
        } else {
            /* reset for bitrate switching */
            setup->attack_handling = 0;

            setup->attdec_filter_mem[0] = 0;
            setup->attdec_filter_mem[1] = 0;

            setup->attdec_detected   = 0;
            setup->attdec_position   = 0;
            setup->attdec_acc_energy = 0;
        }

        bitsTmp = setup->total_bits;
        if (encoder->frame_ms == 2.5) {
            bitsTmp = bitsTmp * 4.0 * (1.0 - 0.4);
        }
        if (encoder->frame_ms == 5) {
            bitsTmp = bitsTmp * 2 - 160;
        }

        if (bitsTmp < 400 + (encoder->fs_idx - 1) * 80) {
            setup->ltpf_enable = 1;
        } else if (bitsTmp < 480 + (encoder->fs_idx - 1) * 80) {
            setup->ltpf_enable = 1;
        } else if (bitsTmp < 560 + (encoder->fs_idx - 1) * 80) {
            setup->ltpf_enable = 1;
        } else if (bitsTmp < 640 + (encoder->fs_idx - 1) * 80) {
            setup->ltpf_enable = 1;
        } else {
            setup->ltpf_enable = 0;
        }
#ifdef ENABLE_HR_MODE
        if (encoder->hrmode) {
            setup->ltpf_enable = 0;
        }

        if (encoder->hrmode && encoder->fs_idx >= 4)
        {
			int real_rate = setup->targetBytes * 8000 / encoder->frame_ms;
			setup->regBits = real_rate / 12500;

			if (encoder->fs_idx == 5)
			{
				if (encoder->frame_ms == 10)
				{
					setup->regBits +=2;
				}
				if (encoder->frame_ms == 2.5)
				{
					setup->regBits -= 6;
				}
			}
			else
			{
				if (encoder->frame_ms == 2.5)
				{
					setup->regBits -= 6;
				}
				else if (encoder->frame_ms == 5)
				{
					setup->regBits += 0;
				}
				if (encoder->frame_ms == 10)
				{
					setup->regBits += 5;
				}
			}
			assert(setup->regBits >= 0);
        }
        else
        {
        	setup->regBits = -1;
        }
#endif
    }

    encoder->bitrate = bitrate;

    return LC3_OK;
}

void update_enc_bandwidth(LC3_Enc* encoder, int bandwidth)
{
    int index = 0;

    if (bandwidth >= encoder->fs_in) {
        encoder->bandwidth = 0;
    } else {
        encoder->bandwidth = bandwidth;
        index              = FS2FS_IDX(bandwidth);
#ifdef ENABLE_HR_MODE
        if (index > 4) {
            index = 5;
        }
#endif
        encoder->bw_ctrl_cutoff_bin = encoder->cutoffBins[index];
    }
}
