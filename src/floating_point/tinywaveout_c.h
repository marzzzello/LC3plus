/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

/*
    DISCLAIMER

    This software module was originally written by Stefan Gewinner and is
    hereby placed in the public domain for all purposes, whether commercial,
    free [as in speech] or educational, etc.

    This code is provided on an "AS IS" basis, WITHOUT WARRANTY OF ANY KIND,
    EITHER EXPRESS OR IMPLIED, AND THE AUTHOR HEREBY DISCLAIMS ALL SUCH
    WARRANTIES, INCLUDING BUT NOT LIMITED TO ANY WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE OR NONINFRINGEMENT.

    The author assumes no liability whatsoever, including infringement of
    any patent or copyright and any other damages, including any general,
    special, incidential or consequential damages that may arise out of
    the use or inability to use this code or any portion thereof (including
    but not limited to loss of data being rendered inaccurate or losses
    sustained by you or third parties or a failure of this code to operate
    with any other programs).

    Copyright (c) 2002 Stefan Gewinner.
    modified to C: 2004 mh

    $Header: /home/cvs/amm/menc/tinywaveout1/include/tinywaveout_c.h,v 1.21
   2014-11-27 10:23:58 jdr Exp $
*/

#ifndef __TINYWAVEOUT_C_H__
#define __TINYWAVEOUT_C_H__

#include <stdio.h>

/*#define SUPPORT_BWF*/

#ifdef SUPPORT_BWF
#include <string.h>
#endif

#if defined(__i386__) || defined(_M_IX86) || defined(_M_X64) ||                \
    defined(__x86_64__) || defined(__arm__) || defined(__aarch64__)
#define __TWO_LE /* _T_iny _W_ave _O_ut _L_ittle _E_ndian */
#endif

#if defined(__POWERPC__)
#define __TWO_BE /* _T_iny _W_ave _O_ut _B_ig _E_ndian */
#endif

#if defined(__sparc__)
#define __TWO_BE /* _T_iny _W_ave _O_ut _B_ig _E_ndian */
#endif

#if !defined(__TWO_LE) && !defined(__TWO_BE)
#error unknown processor
#endif

#define __TWO_SUCCESS (0)
#define __TWO_ERROR (-1)

/*--- local types/structs ----------------------------------*/

#if defined(_MSC_VER)
#pragma pack(push, 1)
#else
#pragma pack(1)
#endif

#ifndef TW_INT64
#if !(defined(WIN32))
#define TWO_INT64 long long
#else
#define TWO_INT64 __int64
#endif
#endif

#ifdef SUPPORT_BWF
typedef struct {
  float loudnessVal;
  float loudnessRange;
  float maxTruePeakLevel;
  float maxMomentaryLoudnes;
  float maxShortTermLoudness;
} LOUDNESSINFO;
#else
typedef void LOUDNESSINFO;
#endif

typedef struct __tinyWaveOutHeader {
  unsigned int riffType;
  unsigned int riffSize;

  unsigned int waveType;
} __tinyWaveOutHeader;

#ifdef SUPPORT_BWF
typedef struct __tinyWaveOutBextChunk {
  unsigned int formatType; /* = 'bext' */
  unsigned int formatSize; /* size info */

  unsigned char description[256];
  unsigned char originator[32];
  unsigned char originatorReference[32];
  unsigned char originatorDate[10]; /* ASCII: <<yyyy:mm:dd>> */
  unsigned char originationTime[8]; /* ASCII: <<hh:mm:ss>> */
  unsigned int timeReferenceLow;
  unsigned int timeReferenceHigh;
  unsigned short version;
  unsigned char UMID[64]; /* Binary Bytes of SMPTE UMID */

  signed short loudnessVal;
  signed short loudnessRange;
  signed short maxTruePeakLevel;
  signed short maxMomentaryLoudnes;
  signed short maxShortTermLoudness;

  unsigned char Reserved[180];

  unsigned char
      codingHistory; /* ASCII: <<History coding>> - undefined length! */
                     /* for variable length, mve this out of this struct */
} __tinyWaveOutBextChunk;
#endif

typedef struct __tinyWaveOutFmtChunk {
  unsigned int formatType;
  unsigned int formatSize;

  unsigned short formatTag;
  unsigned short numChannels;
  unsigned int sampleRate;
  unsigned int bytesPerSecond;
  unsigned short blockAlignment;
  unsigned short bitsPerSample;

  /* wav fmt ext hdr here */
} __tinyWaveOutFmtChunk;

typedef struct __tinyWaveOutDataChunk {
  unsigned int dataType;
  unsigned int dataSize;

} __tinyWaveOutDataChunk;

typedef struct __tinyWaveOutHandle {
  FILE *theFile;
  unsigned int dataSize;
  TWO_INT64 dataSizeLimit;
  unsigned int fmtChunkOffset;
#ifdef SUPPORT_BWF
  unsigned int bextChunkOffset;
#endif
  unsigned int dataChunkOffset;
  unsigned int bps;
  unsigned int clipCount;
} __tinyWaveOutHandle, WAVEFILEOUT;

/*--- local protos --------------------------------------------------*/
static __inline unsigned int BigEndian32(char, char, char, char);
static __inline unsigned int LittleEndian32(unsigned int);
static __inline unsigned int LittleEndian32s(int);
static __inline short LittleEndian16(short);
#ifdef SUPPORT_BWF
static unsigned int EncodeLoudness(float);
#endif
static __inline int __dataSizeChk(WAVEFILEOUT *self, int newbytes);

#if defined(_MSC_VER)
#pragma pack(pop)
#else
#pragma pack()
#endif

#ifdef SUPPORT_BWF
static void setDefaultLoudness(LOUDNESSINFO *x) {
  x->loudnessVal = 1.0f;
  x->loudnessRange = 2.0f;
  x->maxTruePeakLevel = 3.0f;
  x->maxMomentaryLoudnes = 4.0f;
  x->maxShortTermLoudness = 5.0f;
}
#endif

static WAVEFILEOUT *CreateBWF(const char *fileName,
                              const unsigned int sampleRate,
                              const unsigned int numChannels,
                              const unsigned int bps
/* const unsigned int writeWaveExt */
#ifdef SUPPORT_BWF
                              ,
                              const LOUDNESSINFO *hBwfData
#endif
) {
  WAVEFILEOUT *self;
  __tinyWaveOutHeader whdr;
#ifdef SUPPORT_BWF
  __tinyWaveOutBextChunk wbextch;
#endif
  __tinyWaveOutFmtChunk wfch;
  __tinyWaveOutDataChunk wdch;
  unsigned int blockAlignment = 0;
  unsigned int ByteCnt = 0; /* Byte counter for fwrite */

  self = (WAVEFILEOUT *)calloc(1, sizeof(WAVEFILEOUT));
  if (!self)
    goto bail; /* return NULL; */

  if (!fileName)
    goto bail;
  if (sampleRate == 0)
    goto bail;
  if (sampleRate > 192000)
    goto bail;
  if (numChannels == 0)
    goto bail;
  if (numChannels > 48)
    goto bail;
  if (bps != 16 && bps != 24 && bps != 32)
    goto bail;

  self->theFile = fopen(fileName, "wb+");
  if (!self->theFile)
    goto bail;

  /* WAV-Header */
  whdr.riffType = BigEndian32('R', 'I', 'F', 'F');
  whdr.riffSize = LittleEndian32(
      0xffffffff); /* set to maximum, if fseek() doesn't work later */
  whdr.waveType = BigEndian32('W', 'A', 'V', 'E');
  /* write to file */
  ByteCnt = 0;
  ByteCnt += fwrite(&whdr, 1, sizeof(whdr), self->theFile);

#ifdef SUPPORT_BWF
  /* BEXT-Chunk */
  if (hBwfData) {
    memset(&wbextch, 0, sizeof(__tinyWaveOutBextChunk));
    wbextch.formatType = BigEndian32('b', 'e', 'x', 't');
    wbextch.formatSize = LittleEndian32(sizeof(__tinyWaveOutBextChunk) - 8);
    wbextch.version = 0x0002;
    wbextch.loudnessVal = EncodeLoudness(hBwfData->loudnessVal);
    wbextch.loudnessRange = EncodeLoudness(hBwfData->loudnessRange);
    wbextch.maxTruePeakLevel = EncodeLoudness(hBwfData->maxTruePeakLevel);
    wbextch.maxMomentaryLoudnes = EncodeLoudness(hBwfData->maxMomentaryLoudnes);
    wbextch.maxShortTermLoudness =
        EncodeLoudness(hBwfData->maxShortTermLoudness);
    /* t.b.d.:  more values */

    /* write to file */
    self->bextChunkOffset = ByteCnt;
    ByteCnt += fwrite(&wbextch, 1, sizeof(wbextch), self->theFile);
  }
#endif

  /* FMT-Chunk */
  wfch.formatType = BigEndian32('f', 'm', 't', ' ');
  wfch.formatSize = LittleEndian32(16);
  switch (bps) {
  case 16:
  case 24:
    wfch.formatTag = LittleEndian16(0x0001); /* WAVE_FORMAT_PCM */
    break;
  case 32:
    // wfch.formatTag    = LittleEndian16(0x0003);  /* WAVE_FORMAT_IEEE_FLOAT */
    wfch.formatTag = LittleEndian16(0x0001); /* WAVE_FORMAT_PCM */
    break;
  default:
    goto bail;
  }
  self->bps = bps;
  wfch.bitsPerSample = LittleEndian16(bps);
  wfch.numChannels = LittleEndian16(numChannels);
  blockAlignment = numChannels * (bps >> 3);
  wfch.blockAlignment = LittleEndian16(blockAlignment);
  wfch.sampleRate = LittleEndian32(sampleRate);
  wfch.bytesPerSecond = LittleEndian32(sampleRate * blockAlignment);
  /* tbd: wavfmt ext hdr here */
  /* write to file */
  self->fmtChunkOffset = ByteCnt;
  ByteCnt += fwrite(&wfch, 1, sizeof(wfch), self->theFile);

  /* DATA-Chunk */
  self->dataChunkOffset = ByteCnt;
  wdch.dataType = BigEndian32('d', 'a', 't', 'a');
  wdch.dataSize =
      LittleEndian32(0xffffffff - ByteCnt); /* yet unknown. set to maximum */
  /* write to file */
  ByteCnt += fwrite(&wdch, 1, sizeof(wdch), self->theFile);

  self->dataSizeLimit = LittleEndian32(
      0xffffffff - ByteCnt); /* maximum size for data chunk for 4 GB files */
  /* self->dataSizeLimit = LittleEndian32(0x7fffffff - ByteCnt); */ /* maximum
                                                                       size for
                                                                       data
                                                                       chunk for
                                                                       2 GB
                                                                       files */

  self->clipCount = 0;

  return self;

bail:
  free(self);
  return NULL;
}

static WAVEFILEOUT *CreateWav(const char *fileName,
                              const unsigned int sampleRate,
                              const unsigned int numChannels,
                              const unsigned int bps
                              /* const unsigned int writeWaveExt */
) {
  return CreateBWF(fileName, sampleRate, numChannels, bps);
}

#define MAX_PCM16 (+32767)
#define MIN_PCM16 (-32768)
static __inline int CLIP_PCM16(int sample, unsigned int *clipcount) {
  int tmp = sample;

  if (sample >= MAX_PCM16) {
    tmp = MAX_PCM16;
    (*clipcount)++;
  } else {
    if (sample <= MIN_PCM16) {
      tmp = MIN_PCM16;
      (*clipcount)++;
    }
  }

  return tmp;
}

#define MAX_PCM24 (+8388607)
#define MIN_PCM24 (-8388608)
static __inline int CLIP_PCM24(int sample, unsigned int *clipcount) {
  int tmp = sample;

  if (sample >= MAX_PCM24) {
    tmp = MAX_PCM24;
    (*clipcount)++;
  } else {
    if (sample <= MIN_PCM24) {
      tmp = MIN_PCM24;
      (*clipcount)++;
    }
  }

  return tmp;
}

#define MAX_FLOAT32 (+1.0f)
#define MIN_FLOAT32 (-1.0f)
static __inline float CLIP_FLOAT32(float sample, unsigned int *clipcount) {
  float tmp = sample;

  if (sample >= MAX_FLOAT32) {
    tmp = MAX_FLOAT32;
    (*clipcount)++;
  } else {
    if (sample <= MIN_FLOAT32) {
      tmp = MIN_FLOAT32;
      (*clipcount)++;
    }
  }

  return tmp;
}

static int __WriteSample16(WAVEFILEOUT *self, int sample, int scale) {
  size_t cnt;
  short v;

  if ((scale - 16) > 0)
    sample = sample >> (scale - 16);
  else
    sample = sample << (16 - scale);

  v = (short)CLIP_PCM16(sample, &(self->clipCount));
#ifdef __TWO_BE
  v = LittleEndian16(v);
#endif

  cnt = fwrite(&v, sizeof(short), 1, self->theFile);

  if (cnt == 1) {
    self->dataSize += 2;
    return __TWO_SUCCESS;
  }

  return __TWO_ERROR;
}

static int __WriteSample24(WAVEFILEOUT *self, int sample, int scale) {
  size_t cnt;
  int v;

  if ((scale - 24) > 0)
    sample = sample >> (scale - 24);
  else
    sample = sample << (24 - scale);

  v = (int)CLIP_PCM24(sample, &(self->clipCount));
#ifdef __TWO_BE
  v = LittleEndian32s(v);
#endif
  cnt = fwrite(&v, 3, 1, self->theFile);

  if (cnt == 1) {
    self->dataSize += 3;
    return __TWO_SUCCESS;
  }

  return __TWO_ERROR;
}

static int __WriteSample32(WAVEFILEOUT *self, int sample) {
  size_t cnt;
  int v = 0;

  v = sample;
#ifdef __TWO_BE
  v = LittleEndian32s(v);
#endif
  cnt = fwrite(&v, 4, 1, self->theFile);

  if (cnt == 1) {
    self->dataSize += 4;
    return __TWO_SUCCESS;
  }

  return __TWO_ERROR;
}

static int __WriteSampleInt(WAVEFILEOUT *self, int sample, int scale) {
  int err;

  if (!self)
    return __TWO_ERROR;

  switch (self->bps) {

  case 32:
    err = __WriteSample32(self, sample);
    break;

  case 16:
    err = __WriteSample16(self, sample, scale);
    break;

  case 24:
    err = __WriteSample24(self, sample, scale);
    break;

  default:
    err = __TWO_ERROR;
    break;
  }

  return err;
}

/* this function expects values in the 16 bit range +-32767/8 */
/* static int WriteWavShort(
                         WAVEFILEOUT* self,
                         short        sampleBuffer[],
                         unsigned int nSamples
                         )
{
  unsigned long i;
  int err = __TWO_SUCCESS;

  if (!self)         return __TWO_ERROR;
  if (!sampleBuffer) return __TWO_ERROR;
  if (__dataSizeChk(self, nSamples * sizeof(short))) return __TWO_ERROR;

  for (i=0; i< nSamples; i++) {
    if (self->bps == 32)
    {
      err = __WriteSample32(self, sampleBuffer[i] / 32768.0f);
    }
    else
    {
      err = __WriteSampleInt(self, (int)sampleBuffer[i], 16);
    }
    if (err != __TWO_SUCCESS) return err;
  }

  return __TWO_SUCCESS;
}
*/

/* this function expects values in the 24 bit range +-8388607/8 */
static int WriteWavLong(WAVEFILEOUT *self, int sampleBuffer[],
                        unsigned int nSamples) {
  unsigned long i;
  int err = __TWO_SUCCESS;

  if (!self)
    return __TWO_ERROR;
  if (!sampleBuffer)
    return __TWO_ERROR;
  if (__dataSizeChk(self, nSamples * sizeof(int)))
    return __TWO_ERROR;

  for (i = 0; i < nSamples; i++) {
    if (self->bps == 16) {
      err = __WriteSampleInt(self, sampleBuffer[i], 16);
    } else {
      err = __WriteSampleInt(self, sampleBuffer[i], 24);
    }
  }
  if (err != __TWO_SUCCESS)
    return err;

  return __TWO_SUCCESS;
}

/* this function expects normalized values in the range +-1.0 */
#define MAX_FL (+2.0f * 8388608.0f)
#define MIN_FL (-2.0f * 8388608.0f)
#define CLIP_FL(x) (((x) >= MAX_FL) ? MAX_FL : (((x) <= MIN_FL) ? MIN_FL : (x)))
/* static int WriteWavFloat(
                         WAVEFILEOUT* self,
                         float        sampleBuffer[],
                         unsigned int nSamples
                         )
{
  unsigned int i;
  int err = __TWO_SUCCESS;

  if (!self)         return __TWO_ERROR;
  if (!sampleBuffer) return __TWO_ERROR;
  if (__dataSizeChk(self, nSamples * sizeof(float))) return __TWO_ERROR;

  for (i=0; i<nSamples; i++) {
    if(self->bps == 32)
    {
      err = __WriteSample32(self, sampleBuffer[i]);
    }
    else
    {
      float tmp = CLIP_FL(sampleBuffer[i] * 8388608.0f);
      err = __WriteSampleInt(self, (int) tmp, 24);
    }
    if (err != __TWO_SUCCESS) return err;
  }

  return __TWO_SUCCESS;
}
*/

static int CloseWav(WAVEFILEOUT *self) {
  unsigned int riffSize_le = 0;
  unsigned int dataSize_le = 0;

  if (!self)
    return __TWO_ERROR;

  riffSize_le = LittleEndian32(
      self->dataChunkOffset - 8 + 8 +
      self->dataSize); /* sizeof(hdr) - (8 bytes of riff chunk header) + (8
                          bytes data chunk header) + sizeof(raw-pcm-data) */
  dataSize_le = LittleEndian32(self->dataSize);

#ifndef __NO_FSEEK
  /* fseek(self->theFile, 0, SEEK_SET);*/

  /* seek to riffsize */
  fseek(self->theFile, 4, SEEK_SET);
  fwrite(&riffSize_le, sizeof(riffSize_le), 1, self->theFile);
  /* seek to datasize */
  fseek(self->theFile, self->dataChunkOffset + 4, SEEK_SET);
  fwrite(&dataSize_le, sizeof(dataSize_le), 1, self->theFile);

#else
  fclose(self->theFile);
  fopen();
  /* tbd ... */

  /* overwrite the first n bytes w/ wav hdr*/
  fwrite(self->waveHeader, sizeof(__tinyWaveHeader), 1, self->theFile);

#endif

  fclose(self->theFile);
  free(self);

  return __TWO_SUCCESS;
}

#ifdef SUPPORT_BWF
static int CloseBWF(WAVEFILEOUT *self, LOUDNESSINFO bextData) {
  int wordData;

  if (!self)
    return __TWO_ERROR;

  if (self->bextChunkOffset) {
    /* Offset for Loudness Data in bext-chunk: 8: Chunck-Header, 412:prev.Data
     */
    fseek(self->theFile, self->bextChunkOffset + 8 + 412, SEEK_SET);

    wordData = LittleEndian32(EncodeLoudness(bextData.loudnessVal));
    fwrite(&wordData, 2, 1, self->theFile);

    wordData = LittleEndian32(EncodeLoudness(bextData.loudnessRange));
    fwrite(&wordData, 2, 1, self->theFile);

    wordData = LittleEndian32(EncodeLoudness(bextData.maxTruePeakLevel));
    fwrite(&wordData, 2, 1, self->theFile);

    wordData = LittleEndian32(EncodeLoudness(bextData.maxMomentaryLoudnes));
    fwrite(&wordData, 2, 1, self->theFile);

    wordData = LittleEndian32(EncodeLoudness(bextData.maxShortTermLoudness));
    fwrite(&wordData, 2, 1, self->theFile);
  }

  return CloseWav(self);
}
#endif

/*------------- local subs ----------------*/


static __inline unsigned int BigEndian32(char a, char b, char c, char d) {
#ifdef __TWO_LE
  return (unsigned int)d << 24 | (unsigned int)c << 16 | (unsigned int)b << 8 |
         (unsigned int)a;
#else
  return (unsigned int)a << 24 | (unsigned int)b << 16 | (unsigned int)c << 8 |
         (unsigned int)d;
#endif
}

static __inline unsigned int LittleEndian32(unsigned int v) {
#ifdef __TWO_LE
  return v;
#else
  return (v & 0x000000FF) << 24 | (v & 0x0000FF00) << 8 |
         (v & 0x00FF0000) >> 8 | (v & 0xFF000000) >> 24;
#endif
}

/* signed version of the above */
static __inline unsigned int LittleEndian32s(int v) {
#ifdef __TWO_LE
  return v;
#else
  return (v & 0x000000FF) << 24 | (v & 0x0000FF00) << 8 |
         (v & 0x00FF0000) >> 8 | (v & 0xFF000000) >> 24;
#endif
}

static __inline short LittleEndian16(short v) {
#ifdef __TWO_LE
  return v;
#else
  return ((v << 8) & 0xFF00) | ((v >> 8) & 0x00FF);
#endif
}

#ifdef SUPPORT_BWF
static unsigned int EncodeLoudness(float x) {
  int s = (x > 0) - (x < 0);
  return (int)(x * 100.0f + s * 0.5f);
}
#endif

static __inline int __dataSizeChk(WAVEFILEOUT *self, int newbytes) {
  if (!self)
    return __TWO_ERROR;

  if ((((TWO_INT64)self->dataSize) + ((TWO_INT64)newbytes)) >
      self->dataSizeLimit) {
    return __TWO_ERROR;
  }

  return __TWO_SUCCESS;
}

#endif /* __TINYWAVEOUT_C_H__ */
