/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef WAV_IO_H
#define WAV_IO_H

#include <stdio.h>

typedef enum { WAVE_S16 = 16, WAVE_S24 = 24, WAVE_F32 = 32 } WavFormat;

typedef enum {
    WAV_OK           = 0,
    WAV_OPEN_ERROR   = 1,
    WAV_READ_ERROR   = 2,
    WAV_WRITE_ERROR  = 3,
    WAV_SYNTAX_ERROR = 4,
    WAV_FORMAT_ERROR = 5,
    WAV_NULL_ERROR   = 6
} WavError;

typedef struct {
    int       channels;
    int       samplerate;
    WavFormat format;
    int       length;  // in samples
    int       clipped; // total clipped samples
    int       pos;     // internal use
    int       buf_size;
    char*     buf;
    FILE*     file;
} Wav;

// functions
int  wav_reader(Wav* wav, const char* path);
int  wav_read(Wav* wav, float** pcm_buf, int num_samples);
void wav_close(Wav* wav);

#endif
