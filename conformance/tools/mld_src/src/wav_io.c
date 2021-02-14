/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "wav_io.h"
#include "util.h"
#include "shims.h"

#define SAMP16(buf, i, ch, channels) (((uint16_t*)(buf))[(i) * (channels) + (ch)])
#define SAMP32(buf, i, ch, channels) (((uint32_t*)(buf))[(i) * (channels) + (ch)])

static int16_t read_le16(FILE* f)
{
    uint16_t z = 0;
    fread(&z, 1, sizeof(z), f);
    return le16toh(z);
}

static int32_t read_le32(FILE* f)
{
    uint32_t z = 0;
    fread(&z, 1, sizeof(z), f);
    return le32toh(z);
}

static int read_chunk(FILE* f, const char* str)
{
    char chunk[4] = {0};
    fread(chunk, 1, 4, f);
    return !memcmp(chunk, str, 4);
}

int wav_reader(Wav* wav, const char* path)
{
    if (!wav || !path)
        return WAV_NULL_ERROR;
    FILE* f = fopen(path, "rb");
    if (!f)
        return WAV_OPEN_ERROR;
    if (!read_chunk(f, "RIFF")) {
        fclose(f);
        return WAV_SYNTAX_ERROR;
    }
    fseek(f, 4, SEEK_CUR); // skip file size
    if (!read_chunk(f, "WAVE") || !read_chunk(f, "fmt ")) {
        fclose(f);
        return WAV_SYNTAX_ERROR;
    }
    int section_size = read_le32(f);
    int format       = read_le16(f); // 1:pcm 3:float -2:extended
    wav->channels    = read_le16(f);
    wav->samplerate  = read_le32(f);
    fseek(f, 6, SEEK_CUR); // skip byte rate, blockalign
    int sample_size = read_le16(f);
    if (format == -2 && read_le16(f) == 22) {
        sample_size = read_le16(f);
        fseek(f, 4, SEEK_CUR);  // skip channel mask
        format = read_le16(f);  // part of GUID
        fseek(f, 14, SEEK_CUR); // skip rest of GUID
    } else {
        fseek(f, section_size - 16, SEEK_CUR); // skip rest of section
    }
    wav->format = (WavFormat)sample_size;

    bool is_ok  = section_size >= 16 && (format == 1 || format == 3);
    bool is_int = format == 1 && (sample_size == 16 || sample_size == 24);
    bool is_f32 = format == 3 && sample_size == 32;
    if (!(is_ok && (is_int || is_f32))) {
        fclose(f);
        return WAV_FORMAT_ERROR;
    }

    while (!read_chunk(f, "data")) { // skip non data sections
        section_size = read_le32(f);
        if (section_size < 0 || fseek(f, section_size, SEEK_CUR)) {
            fclose(f);
            return WAV_READ_ERROR;
        }
    }
    section_size = read_le32(f);

    wav->length = section_size / (wav->channels * sample_size / 8);
    wav->file   = f;
    return WAV_OK;
}

void wav_close(Wav* wav)
{
    if (!wav || !wav->file)
        return;
    fclose(wav->file);
    free(wav->buf);
    memset(wav, 0, sizeof(*wav));
}

static void* wav_buf(Wav* wav, int size)
{
    if (size > wav->buf_size) {
        wav->buf      = realloc(wav->buf, size);
        wav->buf_size = size;
    }
    return wav->buf;
}

int wav_read(Wav* wav, float** pcm_buf, int num_samples)
{
    const int channels   = wav->channels;
    const int frame_size = channels * wav->format / 8;
    char*     buf        = wav_buf(wav, num_samples * frame_size);

    int n = imin(num_samples, wav->length - wav->pos);
    n     = (int)fread(buf, frame_size, n, wav->file);
    wav->pos += n;

    if (n < num_samples)
        memset(buf + n * frame_size, 0, (num_samples - n) * frame_size);

    for (int ch = 0; ch < channels; ch++) {
        if (wav->format == WAVE_S16) {
            for (int i = 0; i < n; i++) {
                int16_t x      = le16toh(SAMP16(buf, i, ch, channels));
                pcm_buf[ch][i] = x * (1.0 / 32767.0); // not 32768 to keep it conform to matlabs audioread function
            }
        } else if (wav->format == WAVE_S24) {
            for (int i = 0; i < n; i++) {
                int32_t x = 0;
                x |= (unsigned char)buf[3 * i * channels + ch + 0];
                x |= (uint32_t)((unsigned char)buf[3 * i * channels + ch + 1]) << 8;
                x |= (uint32_t)((unsigned char)buf[3 * i * channels + ch + 2]) << 16;
                x |= x & 0x800000 ? 0xff000000 : 0; // extend sign
                pcm_buf[ch][i] = x * (1.0 / 8388607.0); // not 8388608 to keep it conform to matlabs audioread function
            }
        } else if (wav->format == WAVE_F32) {
            for (int i = 0; i < n; i++) {
                uint32_t x     = le32toh(SAMP32(buf, i, ch, channels));
                pcm_buf[ch][i] = *(float*)&x;
            }
        }
    }

    return n;
}
