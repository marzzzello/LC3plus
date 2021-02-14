/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef __TINYWAVEIN_C_H__
#define __TINYWAVEIN_C_H__

/*#define SUPPORT_BWF*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined (__i386__) || defined (_M_IX86) || defined (__x86_64__) || defined (_M_X64) || defined(__arm__)
#define __TWI_LE    /* _T_iny _W_ave _I_n _L_ittle _E_ndian */
#endif

#if defined (__POWERPC__)
#define __TWI_BE    /* _T_iny _W_ave _I_n _B_ig _E_ndian */
#endif

#if ! defined (__TWI_LE) && ! defined (__TWI_BE)
#error unknown processor
#endif


#define __TWI_SUCCESS  (0)
#define __TWI_ERROR    (-1)



#ifdef SUPPORT_BWF
typedef struct {
  float loudnessVal;
  float loudnessRange;
  float maxTruePeakLevel;
  float maxMomentaryLoudnes;
  float maxShortTermLoudness;
} WAVEIN_LOUDNESSINFO;
#endif

typedef struct __tinyWaveInHandle {
  FILE*             theFile;
  fpos_t            dataChunkPos;
  unsigned int      position;
  unsigned int      length;
  unsigned int      bps;
#ifdef SUPPORT_BWF
  WAVEIN_LOUDNESSINFO* loudnessInfo;
#endif
} __tinyWaveInHandle, WAVEFILEIN;


typedef struct
{
  short compressionCode;
  short numberOfChannels ;
  unsigned int  sampleRate ;
  unsigned int  averageBytesPerSecond ;
  short blockAlign ;
  short bitsPerSample ;
  /* short extraFormatBytes ; */
} SWavInfo;

#ifdef SUPPORT_BWF
typedef struct
{
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

  unsigned char codingHistory; /* ASCII: <<History coding>> */
} SBwfWav;
#endif

typedef struct
{
  char chunkID[4] ;
  unsigned int chunkSize ;
  /* long dataOffset ; */ /* never used */
} SChunk;


/* local wrapper, always returns correct endian */
static size_t fread_LE(void *ptr, size_t size, size_t nmemb, FILE *stream);

#ifdef __TWI_BE
static short BigEndian16(short v);
static int BigEndian32(int v);
#endif


/*!
 *  \brief Read header from a WAVEfile. Host endianess is handled accordingly.
 *  \fp filepointer of type FILE*.
 *  \wavinfo SWavInfo struct where the decoded header info is stored into.
 *  \return 0 on success and non-zero on failure.
 *
 */
static WAVEFILEIN* OpenWav(
                           const char* filename,
                           int*        samplerate,
                           int*        channels,
                           int*        samplesInFile,
                           int*        bps
                           )
{
  WAVEFILEIN* self;

  SChunk fmt_chunk, data_chunk;
  int offset;
  int tmpSize ;
  char tmpFormat[4];
  SWavInfo wavinfo = { 0, 0, 0, 0, 0, 0 };


  self = (WAVEFILEIN*)calloc(1, sizeof(WAVEFILEIN));
  if (!self) goto bail; /* return NULL; */

  if (!filename) goto bail;
  if (!samplerate) goto bail;
  if (!channels) goto bail;
  if (!samplesInFile) goto bail;
  if (!bps) goto bail;

  self->theFile = fopen(filename, "rb");
  if (!self->theFile) goto bail;


  /* read RIFF-chunk */
  if (fread(tmpFormat, 1, 4, self->theFile) != 4){
    goto bail;
  }

  if (strncmp("RIFF", tmpFormat, 4)) {
    goto bail;
  }

  /* Read RIFF size. Ignored. */
  fread_LE(&tmpSize, 4, 1, self->theFile);

  /* read WAVE-chunk */
  if (fread(tmpFormat, 1, 4, self->theFile) != 4){
    goto bail;
  }

  if (strncmp("WAVE", tmpFormat, 4)) {
    goto bail;
  }

  /* read format/bext-chunk */
  if (fread(fmt_chunk.chunkID, 1, 4, self->theFile) != 4){
    goto bail;
  }

#ifdef SUPPORT_BWF
  /* test for bext-chunk */
  if (!strncmp("bext", fmt_chunk.chunkID, 4)) {
    /*unsigned int i;*/
    unsigned int bextSize=0;

    if (fread_LE(&bextSize, 1, 4, self->theFile) != 4){
      goto bail;
    }

    self->loudnessInfo = (WAVEIN_LOUDNESSINFO *) calloc(1, sizeof(WAVEIN_LOUDNESSINFO));

    if (bextSize>=602) {       /* minimum size bext-data, w/o 'CodingHistory' */
      int i;
      signed short readBuf=0;
      signed int nulbuf=0;

      /* first skip all descriptive data */
      for(i=0; i<412; i++) {
        if (fread_LE(&nulbuf, 1, 1, self->theFile) != 1){
          goto bail;
        }
        bextSize -=1;
      }
      /* second, read loudness data */
      fread_LE(&readBuf, 2, 1, self->theFile);
      bextSize -=2;
      self->loudnessInfo->loudnessVal = (float)readBuf * 0.01f;

      fread_LE(&readBuf, 2, 1, self->theFile);
      bextSize -=2;
      self->loudnessInfo->loudnessRange = (float)readBuf * 0.01f;

      fread_LE(&readBuf, 2, 1, self->theFile);
      bextSize -=2;
      self->loudnessInfo->maxTruePeakLevel = (float)readBuf * 0.01f;

      fread_LE(&readBuf, 2, 1, self->theFile);
      bextSize -=2;
      self->loudnessInfo->maxMomentaryLoudnes = (float)readBuf * 0.01f;

      fread_LE(&readBuf, 2, 1, self->theFile);
      bextSize -=2;
      self->loudnessInfo->maxShortTermLoudness = (float)readBuf * 0.01f;

      /* skip reserved data */
      for(i=0; i<180; i++) {
        if (fread_LE(&nulbuf, 1, 1, self->theFile) != 1){
          goto bail;
        }
        bextSize -=1;
      }
    }

    /* skip remaining data */
    while (bextSize) {
      int nulbuf;
      if (fread_LE(&nulbuf, 1, 1, self->theFile) != 1){
        goto bail;
      }
      bextSize -=1;
    }


    /* read next chunk header */
    if (fread(fmt_chunk.chunkID, 1, 4, self->theFile) != 4){
      goto bail;
    }
  }
#endif

  /* read a potential JUNK-chunk which is in some WAV files */
  if (!strncmp("JUNK", fmt_chunk.chunkID, 4)) {
    unsigned int junkSize=0;

    if (fread_LE(&junkSize, 1, 4, self->theFile) != 4){
      goto bail;
    }

    /* skip JUNK data */
    while (junkSize) {
      int nulbuf;
      if (fread_LE(&nulbuf, 1, 1, self->theFile) != 1){
        goto bail;
      }
      junkSize -=1;
    }

    /* read next chunk header */
    if (fread(fmt_chunk.chunkID, 1, 4, self->theFile) != 4){
      goto bail;
    }
  }

  /* go on with fmt-chunk */
  if (strncmp("fmt", fmt_chunk.chunkID, 3)) {
    goto bail;
  }

  if (fread_LE(&fmt_chunk.chunkSize, 4, 1, self->theFile) != 1){   /* should be 16 for PCM-format (uncompressed) */
    goto bail;
  }


  /* read  info */
  fread_LE(&(wavinfo.compressionCode),       2, 1, self->theFile);
  fread_LE(&(wavinfo.numberOfChannels),      2, 1, self->theFile);
  fread_LE(&(wavinfo.sampleRate),            4, 1, self->theFile);
  fread_LE(&(wavinfo.averageBytesPerSecond), 4, 1, self->theFile);
  fread_LE(&(wavinfo.blockAlign),            2, 1, self->theFile);
  fread_LE(&(wavinfo.bitsPerSample),         2, 1, self->theFile);

  offset = fmt_chunk.chunkSize - 16;
  if (wavinfo.compressionCode == 0x01) {
    if((wavinfo.bitsPerSample != 16) && (wavinfo.bitsPerSample != 24))
      /* we do only support 16 and 24 bit PCM audio */
      goto bail;
  } else if(wavinfo.compressionCode == 0x03) {
     if(wavinfo.bitsPerSample != 32)
      /* we do only support 32 bit float audio */
      goto bail;
  }

  /* Skip rest of fmt header if any. */
  for (;offset > 0; offset--) {
    fread(&tmpSize, 1, 1, self->theFile);
  }

  do {

    /* Read data chunk ID */
    if (fread(data_chunk.chunkID, 1, 4, self->theFile) != 4){
      goto bail;
    }

    /* Read chunk length. */

    if (fread_LE(&offset, 4, 1, self->theFile) != 1){
      goto bail;
    }

    /* Check for data chunk signature. */
    if (strncmp("data", data_chunk.chunkID, 4) == 0) {
      data_chunk.chunkSize = offset;
      break;
    }

    /* Jump over non data chunk. */
    for (;offset > 0; offset--) {
      fread(&tmpSize, 1, 1, self->theFile);
    }

  } while (!feof(self->theFile));


  /* success so far */
  *samplerate     = wavinfo.sampleRate;
  *channels       = wavinfo.numberOfChannels;
  *samplesInFile  = data_chunk.chunkSize / wavinfo.numberOfChannels;
  *samplesInFile /= ((wavinfo.bitsPerSample + 7) / 8);
  *bps            = wavinfo.bitsPerSample;

  self->position = 0;
  self->bps      = wavinfo.bitsPerSample;
  self->length   = *samplesInFile * wavinfo.numberOfChannels;

  fgetpos(self->theFile, &self->dataChunkPos);

  return self;

 bail:
  free(self);
  return NULL;
}



static int __ReadSample16(
                          WAVEFILEIN* self,
                          int*        sample,
                          int         scale
                          )
{
  size_t cnt;
  short v = 0;


  cnt = fread(&v, 2, 1, self->theFile);

  if (cnt != 1) {
    return __TWI_ERROR;
  }

  self->position += 1;

#ifdef __TWI_BE
  v = BigEndian16(v);
#endif

  if ((scale - 16) > 0)
    *sample = v << (scale - 16);
  else
    *sample = v >> (16 - scale);

  return __TWI_SUCCESS;
}

static int __ReadSample24(
                          WAVEFILEIN* self,
                          int*        sample,
                          int         scale
                          )
{
  size_t cnt;
  int v = 0;


  cnt = fread(&v, 3, 1, self->theFile);

  if (cnt != 1) {
    return __TWI_ERROR;
  }

  self->position += 1;

#ifdef __TWI_BE
  v = BigEndian32(v);
#endif

  if (v >= 0x800000) { v |= 0xff000000; }

  if ((scale - 24) > 0)
    *sample = v << (scale - 24);
  else
    *sample = v >> (24 - scale);

  return __TWI_SUCCESS;
}


static int ReadWavShort(
                        WAVEFILEIN*   self,
                        int           sampleBuffer[],
                        int  nSamplesToRead,
                        int* nSamplesRead
                        )
{
  int  i;
  int err = __TWI_SUCCESS;
  *nSamplesRead = 0;

  if (!sampleBuffer) {
    return __TWI_ERROR;
  }

  /* check if we have enough samples left, if not,
     set nSamplesToRead to number of samples left. */
  if (self->position + nSamplesToRead > self->length) {
    nSamplesToRead = self->length - self->position;
  }

  for (i=0; i< nSamplesToRead; i++) {
    int tmp;
    if (self->bps == 16)
    {
      err = __ReadSample16(self, &tmp, 16);
    } else {
      err = __ReadSample24(self, &tmp, 24);
    }
    if (err != __TWI_SUCCESS) {
      return err;
    }
#if (INPUT_DATA == 16)
    sampleBuffer[i] = (short)tmp;
#else
    sampleBuffer[i] = tmp;
#endif

    *nSamplesRead += 1;
  }

  return __TWI_SUCCESS;
}


static int CloseWavIn(WAVEFILEIN* self)
{
  if (self) {
    if (self->theFile) {
      fclose(self->theFile);
    }
  }
  free(self);

  return __TWI_SUCCESS;
}

/*------------- local subs ----------------*/

static size_t fread_LE(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
#ifdef __TWI_LE
  return fread(ptr, size, nmemb, stream);
#endif
#ifdef __TWI_BE

  unsigned char x[ sizeof(int) ];
  unsigned char* y = (unsigned char*)ptr;
  int i;
  int len;

  len = fread(x, size, nmemb, stream);

  for (i=0; i<size * nmemb; i++) {
    *y++ = x[size * nmemb - i - 1];
  }

  return len;
#endif
}


#ifdef __TWI_BE
static short BigEndian16(short v)
{
  short a = (v & 0x0ff);
  short b = (v & 0x0ff00) >> 8;

  return a <<  8 | b ;
}


static int BigEndian32(int v)
{
  int a = (v & 0x0ff);
  int b = (v & 0x0ff00) >> 8;
  int c = (v & 0x0ff0000) >> 16;
  int d = (v & 0xff000000) >> 24;

  return a << 24 | b << 16 | c <<  8 | d ;
}
#endif

#endif
