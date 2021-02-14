/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include <assert.h>
#include <string.h>     /* for mmove */
#include <stdio.h>   
#include <stdlib.h> 
#include "iisfft.h"
#include "cfft.h"

/*  the fixed length fft functions have been split into sevelral headers to
    have smaller files. to give the compiler more room to optimize the ffts
    can't be in separate compilation units. the header approach seemed to be
    the best compromise. to prevent them being included from anywhere else,
    they are guarded by the INCLUDED_FROM_IISFFT_C macro.
*/
#define INCLUDED_FROM_IISFFT_C
#include "fft_2_9.h"
#include "fft_15_16.h"
#include "fft_32.h"
#include "fft_60_128.h"
#include "fft_240_480.h"
#include "fft_384_768.h"
#include "fft_generic.h"

/* swaps len pairs of values in float vector
 * example: [0 1 2 3] becomes [1 0 3 2] and so forth. */
static void fftf_swapvec(LC3_FLOAT* data, LC3_INT len)
{
    LC3_FLOAT* ptr = data;

    while (ptr < data + 2 * len) {
        LC3_FLOAT tmp = ptr[0];
        ptr[0] = ptr[1];
        ptr[1] = tmp;
        ptr += 2;

    }
}

void LC3_iisfft_apply(Iisfft* handle, LC3_FLOAT* x)
{
  if(handle->sign == -1){
    if (!fft_n(x, handle->length))
        pfaDFT(x, handle->length, handle->scratch, handle->num_factors, handle->factors, handle->scratch2, handle->isPrime);
  }
  else{
      if(!ifft_n(x, handle->length)) {
        fftf_swapvec(x,  handle->length);
        pfaDFT(x, handle->length, handle->scratch, handle->num_factors, handle->factors, handle->scratch2, handle->isPrime);
        fftf_swapvec(x,  handle->length);
      }
  }
}

/* returns 1 if there is no specialized function for length or 1 if a scratch needs to be allocated.
   check the fft_n function */
static LC3_INT need_scratch(LC3_INT n)
{
    return n != 2 && n != 3 && n != 4 && n != 5 && n != 7 && n != 8 && n != 9 && n != 15 && n != 16 && n != 32 &&
           n != 60 && n != 64 && n != 128 && n != 240 && n != 256 && n != 384 && n != 480 && n != 512 && n != 768 &&
           n != 1024;

}

IIS_FFT_ERROR LC3_iisfft_plan(Iisfft* handle, LC3_INT length, LC3_INT sign)
{
    memset(handle, 0, sizeof(Iisfft));
    if (length < 2)
        return IIS_FFT_LENGTH_ERROR;
    handle->length = length;
    handle->sign   = sign;
    if (need_scratch(length)) {
        LC3_INT i, lengthOfPrimeScratch = BORDER_FOR_SECOND_SCRATCH; /* only needed for prime numbers bigger than BORDER_FOR_SECOND_SCRATCH */
        if (!factorize(length, &handle->num_factors, handle->factors, handle->isPrime))
            return IIS_FFT_LENGTH_ERROR;
        handle->scratch = (LC3_FLOAT*)malloc(sizeof(LC3_FLOAT) * 2 * length);
        /* create additional scratch for primeFFT() */
        for (i = 0; i < handle->num_factors; i++) {
            if (handle->isPrime[i] == 1 && handle->factors[i] > lengthOfPrimeScratch) {
                lengthOfPrimeScratch = handle->factors[i];
            }
        }
        if (lengthOfPrimeScratch > BORDER_FOR_SECOND_SCRATCH) {
            handle->scratch2 = (LC3_INT*)malloc(sizeof(LC3_INT) * lengthOfPrimeScratch);
            if (!handle->scratch2)
                return IIS_FFT_MEMORY_ERROR;
        }
        if (!handle->scratch)
            return IIS_FFT_MEMORY_ERROR;
    }

    return IIS_FFT_NO_ERROR;
}

void LC3_iisfft_free(Iisfft* handle)
{
    handle->length = 0;
    if(handle->scratch)
      free(handle->scratch);
    if(handle->scratch2)
      free(handle->scratch2);
}

void LC3_fftf_interleave(const LC3_FLOAT* restrict re, const LC3_FLOAT* restrict im, LC3_FLOAT* restrict out, LC3_INT len)
{
    LC3_INT i = 0;
    for (i = 0; i < len; i++) {
        *out++ = *re++;
        *out++ = *im++;
    }
}

void LC3_fftf_deinterleave(const LC3_FLOAT* restrict in, LC3_FLOAT* restrict re, LC3_FLOAT* restrict im, LC3_INT len)
{
    LC3_INT i = 0;
    for (i = 0; i < len; i++) {
        *re++ = *in++;
        *im++ = *in++;
    }
}

/*  generate sine table needed by LC3_rfft_pre/rfft/post. the table must be freed with iisFree */
LC3_FLOAT* LC3_create_sine_table(LC3_INT len)
{
    LC3_INT i = 0;
    LC3_FLOAT* sine_table = (LC3_FLOAT*)malloc(sizeof(LC3_FLOAT) * (len / 2 + 1));
    if (!sine_table)
        return NULL;

    for (i = 0; i < len / 2 + 1; i++)
        sine_table[i] = (LC3_FLOAT)sin(2.0 * M_PIl * i / len);

    return sine_table;
}

void LC3_rfft_post(const LC3_FLOAT* restrict sine_table, LC3_FLOAT* restrict buf, LC3_INT len)
{
    LC3_FLOAT tmp1, tmp2, tmp3, tmp4, s, c;
    LC3_INT i = 0;

    tmp1 = buf[0] + buf[1];
    buf[1] = buf[0] - buf[1];
    buf[0] = tmp1;

    for (i = 1; i <= (len + 2) / 4; i++) {
        s = sine_table[i];           /* sin(pi*i/(len/2)) */
        c = sine_table[i + len / 4]; /* cos(pi*i/(len/2)) */

        tmp1 = buf[2 * i] - buf[len - 2 * i];
        tmp2 = buf[2 * i + 1] + buf[len - 2 * i + 1];
        tmp3 = s * tmp1 - c * tmp2; /* real part of j*W(k,N)*[T(k) - T'(N-k)] */
        tmp4 = c * tmp1 + s * tmp2; /* imag part of j*W(k,N)*[T(k) - T'(N-k)] */
        tmp1 = buf[2 * i] + buf[len - 2 * i];
        tmp2 = buf[2 * i + 1] - buf[len - 2 * i + 1];

        buf[2 * i] = 0.5f * (tmp1 - tmp3);
        buf[2 * i + 1] = 0.5f * (tmp2 - tmp4);
        buf[len - 2 * i] = 0.5f * (tmp1 + tmp3);
        buf[len - 2 * i + 1] = -0.5f * (tmp2 + tmp4);


    }
}

void LC3_rfft_pre(const LC3_FLOAT* restrict sine_table, LC3_FLOAT* restrict buf, LC3_INT len)
{
    const LC3_FLOAT scale = 1.0f / len;
    LC3_FLOAT tmp1, tmp2, tmp3, tmp4, s, c;
    LC3_INT i = 0;

    tmp1 = buf[0] + buf[1];
    buf[1] = scale * (buf[0] - buf[1]);
    buf[0] = scale * tmp1;

    for (i = 1; i <= (len + 2) / 4; i++) {
        s = sine_table[i];           /* sin(pi*i/(len/2)) */
        c = sine_table[i + len / 4]; /* cos(pi*i/(len/2)) */

        tmp1 = buf[2 * i] - buf[len - 2 * i];
        tmp2 = buf[2 * i + 1] + buf[len - 2 * i + 1];
        tmp3 = s * tmp1 + c * tmp2;  /* real part of j*W(k,N)*[T(k) - T'(N-k)] */
        tmp4 = -c * tmp1 + s * tmp2; /* imag part of j*W(k,N)*[T(k) - T'(N-k)] */
        tmp1 = buf[2 * i] + buf[len - 2 * i];
        tmp2 = buf[2 * i + 1] - buf[len - 2 * i + 1];

        buf[2 * i] = scale * (tmp1 + tmp3);
        buf[2 * i + 1] = -scale * (tmp2 + tmp4);
        buf[len - 2 * i] = scale * (tmp1 - tmp3);
        buf[len - 2 * i + 1] = scale * (tmp2 - tmp4);
    }
}
