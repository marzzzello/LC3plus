/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "lc3.h"
#include "defines.h"
#include "functions.h"
#include <stdio.h>

#include "setup_dec_lc3.h"
#include "setup_enc_lc3.h"

#define RETURN_IF(cond, error) \
    if (cond)                  \
    return (error)

/* ensure api header constants are up to date */
STATIC_ASSERT(LC3_MAX_SAMPLES >= MAX_LEN);
STATIC_ASSERT(LC3_MAX_CHANNELS >= MAX_CHANNELS);
STATIC_ASSERT(LC3_MAX_BYTES >= BYTESBUFSIZE);

/* misc functions ************************************************************/

int lc3_version(void)
{
    return LC3_VERSION;
}

int lc3_channels_supported(int channels)
{
    return channels >= 1 && channels <= MAX_CHANNELS;
}

int lc3_samplerate_supported(int samplerate)
{
    switch (samplerate) {
    case 8000:
        return 1;
    case 16000:
        return 1;
    case 24000:
        return 1;
    case 32000:
        return 1;
    case 44100:
        return 1;
    case 48000:
        return 1;
#ifdef ENABLE_HR_MODE_FLAG
    case 96000:
        return 1;
#endif
    default:
        return 0;
    }
}

static int lc3_plc_mode_supported(LC3_PlcMode plc_mode)
{
    switch ((int)plc_mode)
    {
    case LC3_PLC_STANDARD: return 1;
    default: return 0;
    }
}

static int lc3_frame_size_supported(float frame_ms)
{
    switch ((int)(frame_ms * 10))
    {
    case 25: /* fallthru */
    case 50: /* fallthru */
    case 100: return 1;
    default: return 0;
    }
}

static int null_in_list(void **list, int n)
{
    while (--n >= 0)
        RETURN_IF(list[n] == NULL, 1);
    return 0;
}

/* return pointer to aligned base + base_size, *base_size += size + 4 bytes align */
void *balloc(void *base, size_t *base_size, size_t size)
{
    uintptr_t ptr = ((uintptr_t)base + *base_size + 3) & ~3;
    assert((uintptr_t)base % 4 == 0); /* base must be 4-byte aligned */
    *base_size = (*base_size + size + 3) & ~3;
    return (void *)ptr;
}

/* encoder functions *********************************************************/

LC3_Error lc3_enc_init(LC3_Enc *encoder, int samplerate, int channels)
{
    RETURN_IF(encoder == NULL, LC3_NULL_ERROR);
    RETURN_IF((uintptr_t)encoder % 4 != 0, LC3_ALIGN_ERROR);
    RETURN_IF(!lc3_samplerate_supported(samplerate), LC3_SAMPLERATE_ERROR);
    RETURN_IF(!lc3_channels_supported(channels), LC3_CHANNELS_ERROR);
    return FillEncSetup(encoder, samplerate, channels); /* real bitrate check happens here */
}

int lc3_enc_get_size(int samplerate, int channels)
{
    RETURN_IF(!lc3_samplerate_supported(samplerate), 0);
    RETURN_IF(!lc3_channels_supported(channels), 0);
    return alloc_encoder(NULL, channels);
}

int lc3_enc_get_input_samples(const LC3_Enc *encoder)
{
    RETURN_IF(encoder == NULL, 0);
    return encoder->frame_length;
}

int lc3_enc_get_num_bytes(const LC3_Enc *encoder)
{
    RETURN_IF(encoder == NULL, 0);
    int ch = 0, num_bytes = 0;
    for (ch = 0; ch < encoder->channels; ch ++)
    {
    	num_bytes += encoder->channel_setup[ch]->targetBytes;
    }

    return num_bytes;
}

int lc3_enc_get_real_bitrate(const LC3_Enc *encoder)
{
    int ch = 0, totalBytes = 0;
    RETURN_IF(encoder == NULL, 0);
    RETURN_IF(!encoder->lc3_br_set, LC3_BITRATE_UNSET_ERROR);
    for (ch = 0; ch < encoder->channels; ch++)
    {
        totalBytes += encoder->channel_setup[ch]->targetBytes;
    }
    int bitrate = (totalBytes * 80000)/ encoder->frame_dms;
    if (encoder->fs_in == 44100)
    {
    	int rem = bitrate % 480;
    	bitrate = ((bitrate - rem) / 480)* 441 + (rem * 441) / 480;
    }
    return bitrate;
}

LC3_Error lc3_enc_set_bitrate(LC3_Enc *encoder, int bitrate)
{
    RETURN_IF(encoder == NULL, LC3_NULL_ERROR);
    RETURN_IF(bitrate <= 0, LC3_BITRATE_ERROR);
#ifndef STRIP_HR_MODE_API
#ifdef ENABLE_HR_MODE
    RETURN_IF(encoder->fs_idx == 5 && encoder->hrmode == 0, LC3_HRMODE_ERROR);
#endif
#endif
    return update_enc_bitrate(encoder, bitrate);
}

int lc3_enc_get_delay(const LC3_Enc *encoder)
{
    RETURN_IF(encoder == NULL, 0);
    return encoder->frame_length - 2 * encoder->la_zeroes;
}

LC3_Error lc3_enc_set_frame_ms(LC3_Enc *encoder, float frame_ms)
{
    RETURN_IF(encoder == NULL, LC3_NULL_ERROR);
    RETURN_IF(!lc3_frame_size_supported(frame_ms), LC3_FRAMEMS_ERROR);
    RETURN_IF(encoder->lc3_br_set, LC3_BITRATE_SET_ERROR);
    encoder->frame_dms = (int)(frame_ms * 10);
    encoder->frame_ms = frame_ms;
    set_enc_frame_params(encoder);
    return LC3_OK;
}

#ifndef STRIP_HR_MODE_API
#ifdef ENABLE_HR_MODE
LC3_Error lc3_enc_set_hrmode(LC3_Enc* encoder, int hrmode)
{
    RETURN_IF(encoder == NULL, LC3_NULL_ERROR);
    RETURN_IF(encoder->fs_in < 48000 && hrmode != 0, LC3_SAMPLERATE_ERROR);
    encoder->hrmode = hrmode > 0;
    set_enc_frame_params(encoder);
    return LC3_OK;
}
#endif
#endif

LC3_Error lc3_enc_set_bandwidth(LC3_Enc *encoder, int bandwidth)
{
    RETURN_IF(encoder == NULL, LC3_NULL_ERROR);
#ifdef ENABLE_HR_MODE_FLAG
#ifdef ENABLE_HR_MODE
    RETURN_IF(encoder->hrmode == 1, LC3_HRMODE_BW_ERROR);
#endif
#endif
    LC3_INT effective_fs = encoder->fs_in;
    if (encoder->bandwidth != bandwidth) {
        if (encoder->fs_in > 40000) {
            effective_fs = 40000;
        }
        if ((bandwidth * 2) > effective_fs) {
            return LC3_BW_WARNING;
        }
        else {
            encoder->bandwidth = bandwidth;
            encoder->bw_ctrl_cutoff_bin = (bandwidth * encoder->frame_dms) / 5000;
            encoder->bw_index = MAX(0, (bandwidth / 4000) - 1);
        }
    }
    return LC3_OK;
}

LC3_Error lc3_enc16(LC3_Enc* encoder, int16_t** input_samples, void* output_bytes, int* num_bytes)
{
    return lc3_enc_fl(encoder, (void**)input_samples, 16, output_bytes, num_bytes);
}

LC3_Error lc3_enc24(LC3_Enc* encoder, int32_t** input_samples, void* output_bytes, int* num_bytes)
{
    return lc3_enc_fl(encoder, (void**)input_samples, 24, output_bytes, num_bytes);
}

LC3_Error lc3_enc32(LC3_Enc* encoder, int32_t** input_samples, void* output_bytes, int* num_bytes)
{
    return lc3_enc_fl(encoder, (void**)input_samples, 32, output_bytes, num_bytes);
}


LC3_Error lc3_enc_fl(LC3_Enc* encoder, void** input_samples, int bitdepth, void* output_bytes, int* num_bytes)
{
    RETURN_IF(!encoder || !input_samples || !output_bytes || !num_bytes, LC3_NULL_ERROR);
    RETURN_IF(null_in_list(input_samples, encoder->channels), LC3_NULL_ERROR);
    RETURN_IF(bitdepth != 16 && bitdepth != 24 && bitdepth != 32, LC3_ERROR);
    *num_bytes = Enc_LC3_fl(encoder, input_samples, output_bytes, bitdepth);
    assert(*num_bytes == lc3_enc_get_num_bytes(encoder));
    return LC3_OK;
}

/* decoder functions *********************************************************/

LC3_Error lc3_dec_init(LC3_Dec* decoder, int samplerate, int channels, LC3_PlcMode plc_mode)
{
    RETURN_IF(decoder == NULL, LC3_NULL_ERROR);
    RETURN_IF(!lc3_samplerate_supported(samplerate), LC3_SAMPLERATE_ERROR);
    RETURN_IF(!lc3_channels_supported(channels), LC3_CHANNELS_ERROR);
    RETURN_IF(!lc3_plc_mode_supported(plc_mode), LC3_PLCMODE_ERROR);
    return FillDecSetup(decoder, samplerate, channels, plc_mode);
}

int lc3_dec_get_size(int samplerate, int channels)
{
    RETURN_IF(!lc3_samplerate_supported(samplerate), 0);
    RETURN_IF(!lc3_channels_supported(channels), 0);
    return alloc_decoder(NULL, channels);
}

LC3_Error lc3_dec_set_frame_ms(LC3_Dec* decoder, float frame_ms)
{
    RETURN_IF(decoder == NULL, LC3_NULL_ERROR);
    RETURN_IF(!lc3_frame_size_supported(frame_ms), LC3_FRAMEMS_ERROR);
    RETURN_IF(decoder->plcMeth == 2 && frame_ms != 10, LC3_FRAMEMS_ERROR);
    decoder->frame_ms = frame_ms;
    decoder->frame_dms = (LC3_INT) (frame_ms * 10);
    set_dec_frame_params(decoder);
    return LC3_OK;
}

int lc3_dec_get_output_samples(const LC3_Dec* decoder)
{
    RETURN_IF(decoder == NULL, 0);
    return decoder->frame_length;
}

int lc3_dec_get_delay(const LC3_Dec* decoder)
{
    RETURN_IF(decoder == NULL, 0);
    return decoder->frame_length - 2 * decoder->la_zeroes;
}

LC3_Error lc3_dec_fl(LC3_Dec* decoder, void* input_bytes, int num_bytes, void** output_samples, int bps)
{
    RETURN_IF(!decoder || !input_bytes || !output_samples, LC3_NULL_ERROR);
    RETURN_IF(null_in_list((void**)output_samples, decoder->channels), LC3_NULL_ERROR);
    return Dec_LC3_fl(decoder, input_bytes, num_bytes, output_samples, bps);
}

LC3_Error lc3_dec16(LC3_Dec* decoder, void* input_bytes, int num_bytes, int16_t** output_samples)
{
    return lc3_dec_fl(decoder, input_bytes, num_bytes, (void**)output_samples, 16);
}

LC3_Error lc3_dec24(LC3_Dec* decoder, void* input_bytes, int num_bytes, int32_t** output_samples)
{
    return lc3_dec_fl(decoder, input_bytes, num_bytes, (void**)output_samples, 24);
}

LC3_Error lc3_dec32(LC3_Dec* decoder, void* input_bytes, int num_bytes, int32_t** output_samples)
{
    return lc3_dec_fl(decoder, input_bytes, num_bytes, (void**)output_samples, 32);
}

/* memory functions *********************************************************/

LC3_Error lc3_enc_free_memory(LC3_Enc* encoder)
{
    RETURN_IF(!encoder, LC3_NULL_ERROR);

    lc3_free_encoder_structs(encoder);
    free(encoder);

    return LC3_OK;
}

LC3_Error lc3_dec_free_memory(LC3_Dec* decoder)
{
    RETURN_IF(!decoder, LC3_NULL_ERROR);

    lc3_free_decoder_structs(decoder);
    free(decoder);

    return LC3_OK;
}

LC3_Error lc3_free_encoder_structs(LC3_Enc* encoder)
{
    RETURN_IF(!encoder, LC3_NULL_ERROR);

    int ch = 0;
    for (ch = 0; ch < encoder->channels; ch++) {
        mdct_free(&encoder->channel_setup[ch]->mdctStruct);
        dct2_free(&encoder->channel_setup[ch]->dct2StructSNS);
    }

    return LC3_OK;
}

LC3_Error lc3_free_decoder_structs(LC3_Dec* decoder)
{
    RETURN_IF(!decoder, LC3_NULL_ERROR);

    int ch = 0;
    for (ch = 0; ch < decoder->channels; ch++) {
        dct4_free(&decoder->channel_setup[ch]->dct4structImdct);
    }

    return LC3_OK;
}

#ifndef STRIP_HR_MODE_API_MODE_API
#ifdef ENABLE_HR_MODE
LC3_Error lc3_dec_set_hrmode(LC3_Dec* decoder, int hrmode)
{
    RETURN_IF(decoder == NULL, LC3_NULL_ERROR);
    decoder->hrmode = hrmode > 0;
    set_dec_frame_params(decoder);
    return LC3_OK;
}
#endif
#endif
