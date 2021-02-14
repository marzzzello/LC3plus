/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "functions.h"
#include "setup_enc_lc3.h"

/* if encoder is null only size is reported */
int alloc_encoder(LC3_Enc *encoder, int samplerate, int channels)
{
    int    ch         = 0;
    size_t size       = sizeof(LC3_Enc);
    void * mdct_mem32 = NULL, *stEnc_mdct_mem = NULL;

    for (ch = 0; ch < channels; ch++)
    {
        EncSetup *setup = balloc(encoder, &size, sizeof(EncSetup));
        mdct_mem32      = balloc(encoder, &size, sizeof(*setup->mdct_mem32) * DYN_MAX_MDCT_LEN(samplerate));
        stEnc_mdct_mem  = balloc(encoder, &size, sizeof(*setup->stEnc_mdct_mem) * DYN_MAX_MDCT_LEN(samplerate));
        if (encoder)
        {
            encoder->channel_setup[ch] = setup;
            setup->mdct_mem32          = mdct_mem32;
            setup->stEnc_mdct_mem      = stEnc_mdct_mem;
        }
    }

    return (int)size;
}

LC3_Error FillEncSetup(LC3_Enc *encoder, int samplerate, int channels)
{
    int ch = 0;

    memset(encoder, 0, lc3_enc_get_size(samplerate, channels));
    alloc_encoder(encoder, samplerate, channels);

    encoder->fs                = CODEC_FS(samplerate);
    encoder->fs_in             = samplerate;
    encoder->fs_idx            = FS2FS_IDX(encoder->fs);
    encoder->channels          = channels;
    encoder->frame_dms         = 100;
    encoder->envelope_bits     = 38;
    encoder->global_gain_bits  = 8;
    encoder->noise_fac_bits    = 3;
    encoder->BW_cutoff_bits    = BW_cutoff_bits_all[encoder->fs_idx];
    encoder->r12k8_mem_in_len  = extract_l(L_shr_pos(Mpy_32_16(encoder->fs, 20972), 9));
    encoder->r12k8_mem_out_len = 24;
    encoder->epmr = LC3_EPMR_ZERO;

    for (ch = 0; ch < encoder->channels; ch++)
    {
        encoder->channel_setup[ch]->x_exp      = 15;
        encoder->channel_setup[ch]->resamp_exp = 17;
    }

    set_enc_frame_params(encoder);

    return lc3_enc_set_ep_mode(encoder, LC3_EP_OFF); /* also calls update_enc_bitrate */
}

/* set frame config params */
void set_enc_frame_params(LC3_Enc *encoder)
{
    
    encoder->frame_length = extract_l(L_shr_pos(Mpy_32_16(encoder->fs, 20972), 6)); /* fs * 0.01*2^6 */
    SWITCH (encoder->frame_dms)
    {
    case 25:
        encoder->frame_length       = shr_pos(encoder->frame_length, 2);
        encoder->yLen               = s_min(MAX_BW >> 2, encoder->frame_length);
        encoder->W_fx               = LowDelayShapes_n960_2_5ms[encoder->fs_idx];
        encoder->W_size             = LowDelayShapes_n960_len_2_5ms[encoder->fs_idx];
        encoder->la_zeroes          = LowDelayShapes_n960_la_zeroes_2_5ms[encoder->fs_idx];
        encoder->stEnc_mdct_mem_len = sub(encoder->frame_length, encoder->la_zeroes);
        encoder->bands_number       = bands_number_2_5ms[encoder->fs_idx];
        encoder->bands_offset       = bands_offset_2_5ms[encoder->fs_idx];
        encoder->nSubdivisions      = 2;
        encoder->ltpf_mem_in_len    = LTPF_MEMIN_LEN + 3 * (LEN_12K8 >> 2);
        BREAK;
    case 50:
        encoder->frame_length       = shr_pos(encoder->frame_length, 1);
        encoder->yLen               = s_min(MAX_BW >> 1, encoder->frame_length);
        encoder->W_fx               = LowDelayShapes_n960_5ms[encoder->fs_idx];
        encoder->W_size             = LowDelayShapes_n960_len_5ms[encoder->fs_idx];
        encoder->la_zeroes          = LowDelayShapes_n960_la_zeroes_5ms[encoder->fs_idx];
        encoder->stEnc_mdct_mem_len = sub(encoder->frame_length, encoder->la_zeroes);
        encoder->bands_number       = bands_number_5ms[encoder->fs_idx];
        encoder->bands_offset       = bands_offset_5ms[encoder->fs_idx];
        encoder->nSubdivisions      = 2;
        encoder->ltpf_mem_in_len    = LTPF_MEMIN_LEN + (LEN_12K8 >> 1);
        BREAK;
    case 100:
        encoder->yLen               = s_min(MAX_BW, encoder->frame_length);
        encoder->W_fx               = LowDelayShapes_n960[encoder->fs_idx];
        encoder->W_size             = LowDelayShapes_n960_len[encoder->fs_idx];
        encoder->la_zeroes          = LowDelayShapes_n960_la_zeroes[encoder->fs_idx];
        encoder->stEnc_mdct_mem_len = sub(encoder->frame_length, encoder->la_zeroes);
        encoder->bands_number       = 64;
        encoder->bands_offset       = bands_offset[encoder->fs_idx];
        encoder->nSubdivisions      = 3;
        encoder->ltpf_mem_in_len    = LTPF_MEMIN_LEN;
        BREAK;
    }
}

/* change encoder bitrate */
LC3_Error update_enc_bitrate(LC3_Enc *encoder, int bitrate)
{
    int ch = 0, max_bytes = 0;
    int totalBytes = 0, maxBR = 0, minBR = 0;
    int channel_bytes = 0;

    minBR = MIN_NBYTES * 8 * (10000 / encoder->frame_dms) * encoder->channels   ;
    maxBR = (encoder->fs_in == 44100 ? MAX_NBYTES : MAX_NBYTES_RED) * 8 * (10000 / encoder->frame_dms) * encoder->channels;

    if (bitrate < minBR || bitrate > maxBR)
    {
        return LC3_BITRATE_ERROR;
    }

    encoder->bitrate    = bitrate;
    encoder->lc3_br_set = 1;

    /* move stuff to encoder->channel_setup */

    encoder->combined_channel_coding = 0;
    if (encoder->channels > 1 && encoder->epmode)
    {
        if (encoder->bitrate * encoder->frame_length / (8 * encoder->fs_in) <= 160)
        {
            encoder->combined_channel_coding = 1;
        }
    }

    if (encoder->epmode > 0)
    {
        max_bytes = bitrate * encoder->frame_length / (8 * encoder->fs_in * encoder->channels);
        if (max_bytes < FEC_SLOT_BYTES_MIN || max_bytes > FEC_SLOT_BYTES_MAX)
        {
            return LC3_BITRATE_ERROR;
        }
    }

    if (encoder->combined_channel_coding)
    {
        totalBytes = fec_get_data_size(encoder->epmode, encoder->combined_channel_coding,
                                          bitrate * (Word32)encoder->frame_length / (8 * encoder->fs_in));

        encoder->channel_setup[0]->n_pccw =
            fec_get_n_pccw(bitrate * (Word32)encoder->frame_length / (8 * encoder->fs_in), encoder->epmode,
                           encoder->combined_channel_coding);

        encoder->channel_setup[0]->n_pc = fec_get_n_pc(encoder->epmode, encoder->channel_setup[0]->n_pccw,
                                                       bitrate * (Word32)encoder->frame_length / (8 * encoder->fs_in));
    }
    else
    {
        totalBytes = bitrate * (Word32)encoder->frame_length / (8 * encoder->fs_in);
    }

    for (ch = 0; ch < encoder->channels; ch++)
    {
        EncSetup *setup = encoder->channel_setup[ch];
        channel_bytes = totalBytes / encoder->channels + (ch < (totalBytes % encoder->channels));

        if (encoder->combined_channel_coding)
        {
            setup->targetBytes = channel_bytes;
        }
        else
        {
            setup->targetBytes = fec_get_data_size(encoder->epmode, encoder->combined_channel_coding, channel_bytes);
            setup->n_pccw = fec_get_n_pccw(channel_bytes, encoder->epmode, encoder->combined_channel_coding);
            setup->n_pc = fec_get_n_pc(encoder->epmode, setup->n_pccw, channel_bytes);
        }

        if (encoder->fs_in == 44100)
        {
            max_bytes = MAX_NBYTES;
        }
        else
        {
            max_bytes = MAX_NBYTES_RED;
        }

        if (setup->targetBytes < MIN_NBYTES || setup->targetBytes > max_bytes)
        {
            return LC3_BITRATE_ERROR;
        }

        setup->total_bits = shl(setup->targetBytes, 3);
        setup->targetBitsInit =
            sub(setup->total_bits,
                add(encoder->envelope_bits,
                    add(encoder->global_gain_bits, add(encoder->noise_fac_bits, encoder->BW_cutoff_bits))));
        setup->targetBitsInit = sub(setup->targetBitsInit, sub(17, norm_s(sub(encoder->yLen, 1))));
        if (setup->total_bits > 1280)
        {
            setup->targetBitsInit = sub(setup->targetBitsInit, 1);
        }
        if (setup->total_bits > 2560)
        {
            setup->targetBitsInit = sub(setup->targetBitsInit, 1);
        }

        setup->targetBitsAri = setup->total_bits;

        SWITCH (encoder->frame_dms)
        {
        case 25:
            /* 9830 = 2.4 * 2^12 */
            setup->ltpf_enable =
                sub(extract_l(L_shr(L_mult0(9830, setup->total_bits), 12)), add(560, i_mult(80, encoder->fs_idx))) < 0;
            setup->enable_lpc_weighting = 0;
            BREAK;
        case 50:
            setup->ltpf_enable = sub(sub(i_mult(setup->total_bits, 2), 160), add(560, i_mult(80, encoder->fs_idx))) < 0;
            setup->enable_lpc_weighting = setup->total_bits < 240;
            BREAK;
        case 100:
            setup->enable_lpc_weighting = setup->total_bits < 480;
            setup->ltpf_enable          = sub(setup->total_bits, add(560, i_mult(80, encoder->fs_idx))) < 0;
            BREAK;
        }

        setup->quantizedGainOff =
            -(s_min(115, setup->total_bits / (10 * (encoder->fs_idx + 1))) + 105 + 5 * (encoder->fs_idx + 1));
        if (encoder->frame_dms == 100 && ((encoder->fs_in >= 44100 && setup->targetBytes >= 100) ||
                                          (encoder->fs_in == 32000 && setup->targetBytes >= 81))
#ifdef NONBE_FIX_NO_ATTACK_AT_HIGH_BR
            && setup->targetBytes < 340
#endif
            )
        {
            setup->attack_handling = 1;
        }
        else
        {
            /* reset attack detector for bitrate switching */
            setup->attack_handling      = 0;
            setup->attdec_filter_mem[0] = 0;
            setup->attdec_filter_mem[1] = 0;
            setup->attdec_detected      = 0;
            setup->attdec_position      = 0;
            setup->attdec_acc_energy    = 0;
            setup->attdec_scaling       = 0;
        }
    }

    encoder->bitrate = bitrate;

    return LC3_OK;
}

