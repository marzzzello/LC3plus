/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"
#include "functions.h"
#include "lc3.h"
#include "setup_dec_lc3.h"
#include "setup_enc_lc3.h"

#define RETURN_IF(cond, error)                                                                                         \
    if (cond)                                                                                                          \
    return (error)

#pragma message("PROFILE CONFIG: " PROFILE)
#ifdef SUBSET_NB
#pragma message("- SUBSET_NB")
#endif
#ifdef SUBSET_SQ
#pragma message("- SUBSET_SQ")
#endif
#ifdef SUBSET_HQ
#pragma message("- SUBSET_HQ")
#endif
#ifdef SUBSET_SWB
#pragma message("- SUBSET_SWB")
#endif
#ifdef SUBSET_FB
#pragma message("- SUBSET_FB")
#endif

/* ensure api header constants are up to date */
STATIC_ASSERT(LC3_MAX_SAMPLES >= MAX_LEN);
STATIC_ASSERT(LC3_MAX_CHANNELS >= MAX_CHANNELS);
STATIC_ASSERT(LC3_MAX_BYTES >= BYTESBUFSIZE);
STATIC_ASSERT(LC3_ENC_MAX_SIZE >= ENC_MAX_SIZE);
STATIC_ASSERT(LC3_DEC_MAX_SIZE >= DEC_MAX_SIZE);
STATIC_ASSERT(LC3_ENC_MAX_SCRATCH_SIZE >= SCRATCH_BUF_LEN_ENC_TOT);
STATIC_ASSERT(LC3_DEC_MAX_SCRATCH_SIZE >= SCRATCH_BUF_LEN_DEC_TOT);


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
    switch (samplerate)
    {
#ifdef SUBSET_NB
    case 8000: return 1;
#endif
#ifdef SUBSET_SQ
    case 16000: return 1;
#endif
#ifdef SUBSET_HQ
    case 24000: return 1;
#endif
#ifdef SUBSET_SWB
    case 32000: return 1;
#endif
#ifdef SUBSET_FB
    case 44100: return 1;
    case 48000: return 1;
#endif
    default: return 0;
    }
}

static int lc3_plc_mode_supported(LC3_PlcMode plc_mode)
{
    switch ((int)plc_mode)
    {
    case LC3_PLC_ADVANCED: /* fallthru */
        return 1;
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
    return alloc_encoder(NULL, samplerate, channels);
}

int lc3_enc_get_scratch_size(const LC3_Enc *encoder)
{
    int size = 0;
    RETURN_IF(encoder == NULL, 0);
    size = 14 * MAX(encoder->frame_length, 160) + 64;
    assert(size <= LC3_ENC_MAX_SCRATCH_SIZE);
    return size;
}

int lc3_enc_get_input_samples(const LC3_Enc *encoder)
{
    RETURN_IF(encoder == NULL, 0);
    return encoder->frame_length;
}

int lc3_enc_get_num_bytes(const LC3_Enc *encoder)
{
    RETURN_IF(encoder == NULL, 0);
    return (Word32)encoder->bitrate * encoder->frame_length / (8 * encoder->fs_in);
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
    return update_enc_bitrate(encoder, bitrate);
}

int lc3_enc_get_delay(const LC3_Enc *encoder)
{
    RETURN_IF(encoder == NULL, 0);
    return encoder->frame_length - 2 * encoder->la_zeroes;
}

LC3_Error lc3_enc_set_ep_mode(LC3_Enc *encoder, LC3_EpMode epmode)
{
    RETURN_IF(encoder == NULL, LC3_NULL_ERROR);
    RETURN_IF((unsigned)epmode > LC3_EP_HIGH, LC3_EPMODE_ERROR);
    encoder->epmode = epmode;
    return encoder->lc3_br_set ? update_enc_bitrate(encoder, encoder->bitrate) : LC3_OK;
}

LC3_Error lc3_enc_set_ep_mode_request(LC3_Enc *encoder, LC3_EpModeRequest epmr)
{
    RETURN_IF(encoder == NULL, LC3_NULL_ERROR);
    RETURN_IF((unsigned)epmr > LC3_EPMR_HIGH, LC3_EPMR_ERROR);
    encoder->epmr = epmr;
    return LC3_OK;
}

LC3_Error lc3_enc_set_frame_ms(LC3_Enc *encoder, float frame_ms)
{
    RETURN_IF(encoder == NULL, LC3_NULL_ERROR);
    RETURN_IF(!lc3_frame_size_supported(frame_ms), LC3_FRAMEMS_ERROR);
    RETURN_IF(encoder->lc3_br_set, LC3_BITRATE_SET_ERROR);
    encoder->frame_dms = (int)(frame_ms * 10);
    set_enc_frame_params(encoder);
    return LC3_OK;
}

LC3_Error lc3_enc_set_bandwidth(LC3_Enc *encoder, int bandwidth)
{
    RETURN_IF(encoder == NULL, LC3_NULL_ERROR);
    Word32 effective_fs = encoder->fs_in;
    if (encoder->bandwidth != bandwidth) {
        if (encoder->fs_in > 40000) {
            effective_fs = 40000;
        }
        if ((bandwidth * 2) > effective_fs) {
            return LC3_BW_WARNING;
        }
        else {
            encoder->bandwidth = bandwidth;
            encoder->bw_ctrl_cutoff_bin = (div_l(L_mult0(bandwidth,encoder->frame_dms),(5000>>1)));
            encoder->bw_index = s_max(0,extract_l(div_l(bandwidth,(4000>>1))-1));
        }
    }
    return LC3_OK;
}

static LC3_Error lc3_enc(LC3_Enc *encoder, void **input_samples, int bitdepth, void *output_bytes, int *num_bytes,
                         void *scratch)
{
    RETURN_IF(!encoder || !input_samples || !output_bytes || !num_bytes || !scratch, LC3_NULL_ERROR);
    RETURN_IF(null_in_list(input_samples, encoder->channels), LC3_NULL_ERROR);
    RETURN_IF(bitdepth != 16 && bitdepth != 24, LC3_ERROR);
    RETURN_IF(!encoder->lc3_br_set, LC3_BITRATE_UNSET_ERROR);
    *num_bytes = Enc_LC3(encoder, input_samples, bitdepth, output_bytes, scratch, *num_bytes == -1);
    assert(*num_bytes == lc3_enc_get_num_bytes(encoder));
    return LC3_OK;
}

LC3_Error lc3_enc16(LC3_Enc *encoder, int16_t **input_samples, void *output_bytes, int *num_bytes, void *scratch)
{
    return lc3_enc(encoder, (void **)input_samples, 16, output_bytes, num_bytes, scratch);
}

LC3_Error lc3_enc24(LC3_Enc *encoder, int32_t **input_samples, void *output_bytes, int *num_bytes, void *scratch)
{
    return lc3_enc(encoder, (void **)input_samples, 24, output_bytes, num_bytes, scratch);
}

/* decoder functions *********************************************************/

LC3_Error lc3_dec_init(LC3_Dec *decoder, int samplerate, int channels, LC3_PlcMode plc_mode)
{
    RETURN_IF(decoder == NULL, LC3_NULL_ERROR);
    RETURN_IF(!lc3_samplerate_supported(samplerate), LC3_SAMPLERATE_ERROR);
    RETURN_IF(!lc3_channels_supported(channels), LC3_CHANNELS_ERROR);
    RETURN_IF(!lc3_plc_mode_supported(plc_mode), LC3_PLCMODE_ERROR);
    return FillDecSetup(decoder, samplerate, channels, plc_mode);
}

int lc3_dec_get_size(int samplerate, int channels, LC3_PlcMode plc_mode)
{
    RETURN_IF(!lc3_samplerate_supported(samplerate), 0);
    RETURN_IF(!lc3_channels_supported(channels), 0);
    RETURN_IF(!lc3_plc_mode_supported(plc_mode), LC3_PLCMODE_ERROR);
    return alloc_decoder(NULL, samplerate, channels, plc_mode);
}

int lc3_dec_get_scratch_size(const LC3_Dec *decoder)
{
    int size = 0;
    RETURN_IF(decoder == NULL, 0);
    size = 12 * DYN_MAX_LEN(decoder->fs) + 752;
    if (decoder->plcMeth != LC3_PLC_STANDARD)
        size += 2 * MAX_LGW + 8 * DYN_MAX_LPROT(decoder->fs) + 8 * DYN_MAX_LEN(decoder->fs);
    assert(size <= LC3_DEC_MAX_SCRATCH_SIZE);
    return size;
}

LC3_Error lc3_dec_set_ep_enabled(LC3_Dec *decoder, int ep_enabled)
{
    RETURN_IF(decoder == NULL, LC3_NULL_ERROR);
    decoder->ep_enabled = ep_enabled != 0;
    decoder->epmr       = LC3_EPMR_ZERO;
    return LC3_OK;
}

int lc3_dec_get_error_report(const LC3_Dec *decoder)
{
    RETURN_IF(decoder == NULL, 0);
    return decoder->error_report;
}

LC3_EpModeRequest lc3_dec_get_ep_mode_request(const LC3_Dec *decoder)
{
    RETURN_IF(decoder == NULL, LC3_EPMR_ZERO);
    return (LC3_EpModeRequest)decoder->epmr;
}

LC3_Error lc3_dec_set_frame_ms(LC3_Dec *decoder, float frame_ms)
{
    RETURN_IF(decoder == NULL, LC3_NULL_ERROR);
    RETURN_IF(!lc3_frame_size_supported(frame_ms), LC3_FRAMEMS_ERROR);
    RETURN_IF(decoder->plcMeth == 2 && frame_ms != 10, LC3_FRAMEMS_ERROR);

    decoder->frame_dms = (int)(frame_ms * 10);
    set_dec_frame_params(decoder);
    return LC3_OK;
}

int lc3_dec_get_output_samples(const LC3_Dec *decoder)
{
    RETURN_IF(decoder == NULL, 0);
    return decoder->frame_length;
}

int lc3_dec_get_delay(const LC3_Dec *decoder)
{
    RETURN_IF(decoder == NULL, 0);
    return decoder->frame_length - 2 * decoder->la_zeroes;
}

static LC3_Error lc3_dec(LC3_Dec *decoder, void *input_bytes, int num_bytes, void **output_samples, int bitdepth,
                         void *scratch, int bfi_ext)
{
    RETURN_IF(!decoder || !input_bytes || !output_samples || !scratch, LC3_NULL_ERROR);
    RETURN_IF(null_in_list(output_samples, decoder->channels), LC3_NULL_ERROR);
    RETURN_IF(bitdepth != 16 && bitdepth != 24, LC3_ERROR);
    return Dec_LC3(decoder, input_bytes, num_bytes, output_samples, bitdepth, scratch, bfi_ext);
}

LC3_Error lc3_dec16(LC3_Dec *decoder, void *input_bytes, int num_bytes, int16_t **output_samples, void *scratch, int bfi_ext)
{
    return lc3_dec(decoder, input_bytes, num_bytes, (void **)output_samples, 16, scratch, bfi_ext);
}

LC3_Error lc3_dec24(LC3_Dec *decoder, void *input_bytes, int num_bytes, int32_t **output_samples, void *scratch, int bfi_ext)
{
    return lc3_dec(decoder, input_bytes, num_bytes, (void **)output_samples, 24, scratch, bfi_ext);
}
