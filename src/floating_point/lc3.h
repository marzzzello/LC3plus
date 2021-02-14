/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

/*! \file lc3.h
 *  This header provides the API for LC3plus.
 *
 *  This library is targeting devices with extreme memory limitations, so memory management
 *  must be handeled by the user. This includes allocating memory for the structs. The structs are persistent 
 *  between function calls.
 *
 *  The amount of memory needed for various configurations can be obtained from the lc3_*_get_size
 *  function. The LC3_*_MAX_SIZE macro can be used for all configurations.
 *
 *  Depending on the build configuration some functions might not be available.
 */

#ifndef LC3_H
#define LC3_H

#ifndef _MSC_VER
#include <stdint.h>
#else
typedef unsigned char uint8_t;
typedef __int16       int16_t;
typedef __int32       int32_t;
#endif

/*! Construct version number from major/minor/micro values. */
#define LC3_VERSION_INT(major, minor, micro) (((major) << 16) | ((minor) << 8) | (micro))

/*! Version number to ensure header and binary are matching. */
#define LC3_VERSION LC3_VERSION_INT(1, 4, 2)

/*! Maximum number of supported channels. The actual binary might support
 *  less, use lc3_channels_supported() to check. */
#define LC3_MAX_CHANNELS 16

/*! Maximum number of samples per channel that can be stored in one LC3 frame.
 */
#define LC3_MAX_SAMPLES 960

/*! Maximum number of bytes of one LC3 frame. */
#define LC3_MAX_BYTES 6960

/*! Error codes returned by functions. */
typedef enum {
    LC3_OK                  = 0,  /*!< No error occurred */
    LC3_ERROR               = 1,  /*!< Function call failed */
    LC3_DECODE_ERROR        = 2,  /*!< Frame failed to decode and was concealed */
    LC3_NULL_ERROR          = 3,  /*!< Pointer argument is null */
    LC3_SAMPLERATE_ERROR    = 4,  /*!< Invalid samplerate value */
    LC3_CHANNELS_ERROR      = 5,  /*!< Invalid channels value */
    LC3_BITRATE_ERROR       = 6,  /*!< Invalid bitrate value */
    LC3_NUMBYTES_ERROR      = 7,  /*!< Invalid num_bytes value */
    LC3_EPMODE_ERROR        = 8,  /*!< Invalid ep_mode value */
    LC3_FRAMEMS_ERROR       = 9,  /*!< Invalid frame ms value */
    LC3_ALIGN_ERROR         = 10, /*!< Unaligned pointer */
    LC3_HRMODE_ERROR        = 11, /*!< Invalid usage of hrmode, sampling rate and frame size */
    LC3_BITRATE_UNSET_ERROR = 12, /*!< Function called before bitrate has been set */
    LC3_BITRATE_SET_ERROR   = 13, /*!< Function called after bitrate has been set */
    LC3_HRMODE_BW_ERROR     = 14, /*!< High quality mode and bandwidth switching must not be used together */
    LC3_PLCMODE_ERROR       = 15, /*!< Invalid plc_method value */
        
    /* START WARNING */
    LC3_WARNING             = 16,
    LC3_BW_WARNING          = 17  /*!< Invalid bandwidth cutoff frequency */
} LC3_Error;

/*! Decoder packet loss concealment mode */
typedef enum
{
    LC3_PLC_STANDARD = 0, /*!< Less complex than advanced method */
    LC3_PLC_ADVANCED = 1  /*!< Enhanced concealment method */
} LC3_PlcMode;

typedef struct LC3_Enc LC3_Enc; /*!< Opaque encoder struct. */
typedef struct LC3_Dec LC3_Dec; /*!< Opaque decoder struct. */

/*! \addtogroup Misc
 *  \{ */

/*! Return library version number. It should match LC3_VERSION. */
int lc3_version(void);

/*! Tests if the library supports number of channels.
 *
 *  \param[in]  channels    Number of channels.
 *  \return                 1 for true, 0 for false.
 */
int lc3_channels_supported(int channels);

/*! Tests if the library supports a sampling rate.
 *
 *  \param[in]  samplerate  Sampling rate
 *  \return                 1 for true, 0 for false
 */
int lc3_samplerate_supported(int samplerate);

/*! \}
 *  \addtogroup Encoder
 *  \{ */

/*!
 *  Initialize LC3 encoder.
 *
 *  This function is used to fill a user-allocated encoder struct. This is
 *  typically called once for a samplerate / channel configuration. The bitrate
 *  can be changed later with lc3_enc_set_bitrate().
 *
 *  Recommended bitrates for input sampling rates:
 *     8 kHz: 24 kbps
 *    16 kHz: 32 kbps
 *    24 kHz: 48 kbps
 *    32 kHz: 64 kbps
 *    44.1/48 kHz: 80 kbps (voice), 128 kbps (music)
 *    96 kHz: 156 kbps
 *
 *  \param[out] encoder     Pointer to allocated encoder memory. It must have a
 * size provided by lc3_enc_get_size() for matching samplerate / channels
 *                          configuration or LC3_ENC_MAX_SIZE.
 *  \param[in]  channels    Number of channels.
 *  \param[in]  samplerate  Input sampling rate. Allowed sampling rates are:
 *                          8000, 16000, 24000, 32000, 44100, 48000, 96000
 *  \return                 LC3_OK on success or appropriate error code.
 */
LC3_Error lc3_enc_init(LC3_Enc* encoder, int samplerate, int channels);

/*!
 *  Encode LC3 frame.
 *
 *  Each call consumes a fixed number of samples. The number of input samples
 *  can be obtained from lc3_enc_get_input_samples().
 *
 *  \param[in]  encoder         Encoder handle initialized by lc3_enc_init().
 *  \param[in]  input_samples   Input samples. The left channel is stored in
 *                              input_samples[0], the right channel in input_samples[1]. The input is not
 *                              changed by the encoder.
 *  \param[in]  bps             Number of bits per sample.
 *  \param[out] output_bytes    Output buffer. It must have a at least
 *                              lc3_enc_get_num_bytes() or at most LC3_MAX_BYTES.
 *  \param[out] num_bytes       Number of bytes written to output_bytes.
 *  \return                     LC3_OK on success or appropriate error code.
 */

LC3_Error lc3_enc_fl(LC3_Enc* encoder, void** input_samples, int bitdepth, void* output_bytes, int* num_bytes);

/*! Encode LC3 frame with 16 bit input. See lc3_enc_fl(). */
LC3_Error lc3_enc16(LC3_Enc* encoder, int16_t** input_samples, void* output_bytes, int* num_bytes);

/*! Encode LC3 frame with 24 bit input. See lc3_enc16(). */
LC3_Error lc3_enc24(LC3_Enc* encoder, int32_t** input_samples, void* output_bytes, int* num_bytes);

/*! Encode LC3 frame with 32 bit input. See lc3_enc16(). */
LC3_Error lc3_enc32(LC3_Enc* encoder, int32_t** input_samples, void* output_bytes, int* num_bytes);

/*! Get the size of the LC3 encoder struct for a samplerate / channel
 * configuration. If memory is not restricted LC3_ENC_MAX_SIZE can be used for
 * all configurations.
 *
 *  \param[in]  samplerate  Sampling rate.
 *  \param[in]  channels    Number of channels.
 *  \return                 Size in bytes or 0 on error.
 */
int lc3_enc_get_size(int samplerate, int channels);

/*! Get number of samples per channel expected by lc3_enc16() or lc3_enc24() or lc3_enc32().
 *
 *  \param[in]  encoder     Encoder handle.
 *  \return                 Number of samples or 0 on error.
 */
int lc3_enc_get_input_samples(const LC3_Enc* encoder);

/*! Get real internal bitrate of the encoder. It might differ from the requested
 *  bitrate due to 44.1 kHz input.
 *
 *  \param[in]  encoder     Encoder handle.
 *  \return                 Bitrate in bits per second or 0 on error.
 */
int lc3_enc_get_real_bitrate(const LC3_Enc* encoder);

/*! Set encoder bitrate for all channels.
 *  This function must be called at least once before encoding the first frame, but
 *  after other configuration functions such as lc3_enc_set_frame_ms().
 *
 *  Recommended bitrates for input sampling rates with 10 ms framing:
 *  kHz     | kbps
 *  --------|-----
 *  8       | 24
 *  16      | 32
 *  24      | 48
 *  32      | 64
 *  44.1/48 | 80(voice) 128(music)
 *  96      | 128 
 *
 *  \param[in]  encoder     Encoder handle.
 *  \param[in]  bitrate     Bitrate in bits per second.
 *  \return                 LC3_OK on success or appropriate error code.
 */
LC3_Error lc3_enc_set_bitrate(LC3_Enc* encoder, int bitrate);

/*! Get the maximum number of bytes produced by lc3_enc16() or lc3_enc24() or lc3_enc32() for
 * the current bitrate. It should be equal to the num_bytes output of
 * lc3_enc16() / lc3_enc24() / lc3_enc32().
 *
 *  \param[in]  encoder     Encoder handle.
 *  \return                 Size in bytes or 0 on error.
 */
int lc3_enc_get_num_bytes(const LC3_Enc* encoder);

/*! Set the frame length for LC3 encoder. Allowed values are 10 (default), 5
 * ms and 2.5 ms. The decoder must be configured with lc3_dec_set_frame_ms() with the same
 * value.
 *
 *  \param[in]  encoder     Encoder handle.
 *  \param[in]  frame_ms    Frame length in ms.
 *  \return                 LC3_OK on success or appropriate error code.
 */
LC3_Error lc3_enc_set_frame_ms(LC3_Enc* encoder, float frame_ms);

/*! Set the high resolution mode for LC3 encoder. This mode is mandatory for 96 kHz input and can
 *  also be used for 48 kHz input. Encoder and decoder 
 *  must have the same high resolution mode active.
 *
 *  \param[in]  encoder     Encoder handle.
 *  \param[in]  hrmode      High resolution mode either 1 or 0.
 *  \return                 LC3_OK on success or appropriate error code.
 */
LC3_Error lc3_enc_set_hrmode(LC3_Enc* encoder, int hrmode);

/*! Set encoder bandwidth to a different value. All frequency bins above the cutoff
 *  frequency are cut off. Allowed frequencies are: 4 kHz, 8 kHz, 12 kHz, 16 kHz and 24 kHz.
 *
 *  \param[in]  encoder     Encoder handle.
 *  \param[in]  bandwidth   Cutoff Frequency in Hz
 *  \return                 LC3_OK on success or appropriate error code.
 */
LC3_Error lc3_enc_set_bandwidth(LC3_Enc* encoder, int bandwidth);

/*! Get the encoder delay in number of samples.
 *
 *  \param[in]  encoder     Encoder handle.
 *  \return                 Encoder in samples or 0 on error.
 */
int lc3_enc_get_delay(const LC3_Enc* encoder);

/*! Free memory allocated within lc3 encoder struct.
 *
 *  \param[in]  encoder     Encoder handle.
 *  \return                 LC3_OK on success or appropriate error code.
 */
LC3_Error lc3_enc_free_memory(LC3_Enc* encoder);

/*! Internal function called by lc3_enc_free_memory.
 *
 *  \param[in]  encoder     Encoder handle.
 *  \return                 LC3_OK on success or appropriate error code.
 */
LC3_Error lc3_free_encoder_structs(LC3_Enc* encoder);

/*! \}
 *  \addtogroup Decoder
 *  \{ */

/*!
 *  Initialize LC3 decoder.
 *
 *  This function is used to fill a user-allocated decoder struct. This is
 * typically called once for a samplerate / channel configuration.
 *
 *  The samplerate and channel arguments must have the same values that were
 * used for encoding. LC3 does not provide a signalling scheme, transporting
 * these values is the responsibility of the application.
 *
 *  \param[out] decoder         Pointer to decoder memory. It must have as size
 *                              of least lc3_dec_get_size() or at most LC3_DEC_MAX_SIZE.
 *  \param[in] samplerate       Bitstream sampling rate. \param[in]  channels Bitstream
 *                              number of channels.
 *
 *  \return                     LC3_OK on success or appropriate error code.
 */
LC3_Error lc3_dec_init(LC3_Dec* decoder, int samplerate, int channels, LC3_PlcMode plc_mode);

/*!
 *  Decode LC3 frame with bps bit output.
 *
 *  Each call decodes a fixed number of samples. Use
 * lc3_dec_get_output_samples() to obtain this number. When the input is
 * corrupted and can not be decoded, LC3_DECODE_ERROR is returned and packet
 * loss concealment is applied, so the output is still usable.
 *
 *  \param[in]  decoder         Decoder initialized by lc3_dec_init().
 *  \param[in]  input_bytes     Input bytes.
 *  \param[in]  bps             Bits per audio sample for correct scaling in the end.
 *  \param[in]  num_bytes       Number of valid bytes in input_bytes. To signal a lost frame and
 *                              generate concealment output this value must be set to 0.
 *  \param[out] output_samples  Array of pointers to output channel buffers. Each channel buffer should provide
 *                              enough space to hold at most LC3_MAX_SAMPLES. The left channel is stored in
 *                              output_samples[0], the right channel in output_samples[1].
 *  \return                     Returns LC3_OK on success or appropriate error code. Note
 *                              there is a special case for LC3_DECODE_ERROR where the output is still valid.
 */

LC3_Error lc3_dec_fl(LC3_Dec* decoder, void* input_bytes, int num_bytes, void** output_samples, int bps);

/*! Decode LC3 frame with 16 bit output. See lc3_dec_fl(). */
LC3_Error lc3_dec16(LC3_Dec* decoder, void* input_bytes, int num_bytes, int16_t** output_samples);

/*! Decode LC3 frame with 24 bit output. See lc3_dec_fl(). */
LC3_Error lc3_dec24(LC3_Dec* decoder, void* input_bytes, int num_bytes, int32_t** output_samples);

/*! Decode LC3 frame with 32 bit output. See lc3_dec_fl(). */
LC3_Error lc3_dec32(LC3_Dec* decoder, void* input_bytes, int num_bytes, int32_t** output_samples);

/*! Get the size of the LC3 decoder struct for a samplerate / channel
 * configuration. If memory is not restricted LC3_DEC_MAX_SIZE can be used for
 * all configurations.
 *
 *  \param[in]  channels    Number of channels.
 *  \param[in]  samplerate  Sampling rate.
 *  \return                 Size in bytes or 0 on error.
 */
int lc3_dec_get_size(int samplerate, int channels);

/*! Set the frame length for LC3 decoder. Allowed values are 10 (default) 5
 * ms and 2.5 ms. This only works if the encoder was configured with the same vale.
 *
 *  \param[in]  decoder     Decoder handle.
 *  \param[in]  frame_ms    Frame length in ms.
 *  \return                 LC3_OK on success or appropriate error code.
 */
LC3_Error lc3_dec_set_frame_ms(LC3_Dec* decoder, float frame_ms);

/*! Get the number of samples per channel produced by lc3_dec16() or lc3_dec24 or lc3_dec32().
 *
 *  \param[in]  decoder     Decoder handle.
 *  \return                 Number of samples or 0 on error.
 */
int lc3_dec_get_output_samples(const LC3_Dec* decoder);

/*! Get the decoder delay in number of samples.
 *
 *  \param[in]  decoder     Decoder handle.
 *  \return                 Delay in samples or 0 on error.
 */
int lc3_dec_get_delay(const LC3_Dec* decoder);

/*! Set the high resolution for LC3 decoder. This mode is mandatory for 96 kHz input and can
 *  also be used for 48 kHz input. Encoder and decoder 
 *  must have the same high resolution mode active.
 *
 *  \param[in]  encoder     Decoder handle.
 *  \param[in]  hrmode      High resolution mode either 1 or 0.
 *  \return                 LC3_OK on success or appropriate error code.
 */
LC3_Error lc3_dec_set_hrmode(LC3_Dec* decoder, int hrmode);

/*! Free memory allocated within lc3 decoder struct.
 *
 *  \param[in]  decoder     Decoder handle.
 *  \return                 LC3_OK on success or appropriate error code.
 */
LC3_Error lc3_dec_free_memory(LC3_Dec* decoder);

/*! Internal function called by lc3_dec_free_memory.
 *
 *  \param[in]  decoder     Decoder handle.
 *  \return                 LC3_OK on success or appropriate error code.
 */
LC3_Error lc3_free_decoder_structs(LC3_Dec* decoder);

/*! \} */
#endif /* LC3 */
