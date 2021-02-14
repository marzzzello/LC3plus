/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include <assert.h>
#include <stddef.h>
#include "iis_fft.h"
#include <stdio.h>   
#include <stdlib.h>
#include <string.h>
#include "../structs.h"

/**************************************************************************************************/

/* AFFT uses two fft implementations
 * cfft is used for lengths of power of two, >= 256.
 * iisfft is used for everything else. it is optimized for certain lengths. for a list of
   fast lengths, check the fft_n function.
*/

#include <math.h>
#include "cfft.h"
#include "iisfft.h"


#define FFT_COMPLEX 1
#define FFT_REAL 2

typedef struct T_IIS_FFT
{
    IIS_FFT_DIR sign;
    LC3_INT len;
    LC3_FLOAT* buffer;
    LC3_FLOAT* sine_table;
    Iisfft iisfft;
    Cfft cfft;
} IIS_FFT;


static IIS_FFT_ERROR create(HANDLE_IIS_FFT* handle, LC3_INT type, LC3_INT len, IIS_FFT_DIR sign)
{
    IIS_FFT_ERROR err = IIS_FFT_MEMORY_ERROR;
    HANDLE_IIS_FFT h = NULL;
    /* for real transforms the actual performed fft is half length */
    LC3_INT trlen = (type == FFT_COMPLEX) ? len : len / 2;

    /* check argument sanity */
    if (sign != IIS_FFT_FWD && sign != IIS_FFT_BWD)
        return IIS_FFT_INTERNAL_ERROR;

    
    if (len < 2 || (type == FFT_REAL && len % 4 != 0))
        return IIS_FFT_LENGTH_ERROR;

    h = (HANDLE_IIS_FFT)calloc(1, sizeof(IIS_FFT));
    if (!h)
        return IIS_FFT_MEMORY_ERROR;

    h->len = len;
    h->sign = sign;

    /* temp buffer is not needed for cfft, which is only used for length of n^2 and >= 256 */
    if (!(trlen >= 256 && CFFT_PLAN_SUPPORT(trlen) && type == FFT_COMPLEX)) {
        h->buffer = (LC3_FLOAT*)malloc(sizeof(LC3_FLOAT) * trlen * 2);
        if (!h->buffer)
            goto handle_error1;
    }

    /* create sine lookup table for real ffts */
    if (type == FFT_REAL) {
        h->sine_table = LC3_create_sine_table(len);
        if (!h->sine_table)
            goto handle_error1;
    }

    /* use cfft for legth of power two larger than 256. for length below iisfft is faster */
    if (trlen >= 256 && CFFT_PLAN_SUPPORT(trlen)) {
        LC3_INT s = (type == FFT_REAL) ? IIS_FFT_FWD : sign;
        err = LC3_cfft_plan(&h->cfft, trlen, s) ? IIS_FFT_NO_ERROR : IIS_FFT_INTERNAL_ERROR;
    } else {
        LC3_INT s = (type == FFT_REAL) ? IIS_FFT_FWD : sign;
        err = LC3_iisfft_plan(&h->iisfft, trlen, s);
    }
        if (err != IIS_FFT_NO_ERROR)
            goto handle_error2;

    *handle = h;
    return IIS_FFT_NO_ERROR;

handle_error2:
    free(h->buffer);
handle_error1:
    free(h);
    
    return err;
}

static IIS_FFT_ERROR destroy(HANDLE_IIS_FFT* handle)
{
    if (handle && *handle) {
        LC3_iisfft_free(&(*handle)->iisfft);
        LC3_cfft_free(&(*handle)->cfft);
        if((*handle)->sine_table)
          free((*handle)->sine_table);
        if((*handle)->buffer)
          free((*handle)->buffer);
        free(*handle);
        *handle = NULL;
    }
    return IIS_FFT_NO_ERROR;
}

IIS_FFT_ERROR LC3_IIS_CFFT_Create(HANDLE_IIS_FFT* handle, LC3_INT len, IIS_FFT_DIR sign)
{
    return create(handle, FFT_COMPLEX, len, sign);
}

IIS_FFT_ERROR LC3_IIS_RFFT_Create(HANDLE_IIS_FFT* handle, LC3_INT len, IIS_FFT_DIR sign)
{
    return create(handle, FFT_REAL, len, sign);
}

IIS_FFT_ERROR LC3_IIS_xFFT_Destroy(HANDLE_IIS_FFT* handle) { return destroy(handle); }

IIS_FFT_ERROR LC3_IIS_CFFT_Destroy(HANDLE_IIS_FFT* handle) { return destroy(handle); }

IIS_FFT_ERROR LC3_IIS_RFFT_Destroy(HANDLE_IIS_FFT* handle) { return destroy(handle); }

IIS_FFT_ERROR LC3_IIS_FFT_Apply_CFFT(HANDLE_IIS_FFT handle, const Complex* input, Complex* output)
{
    if (!handle)
        return IIS_FFT_INTERNAL_ERROR;

    /* check for inplace operation */
    memmove(output, input, sizeof(*input) * handle->len);
    LC3_FLOAT* dummy = (LC3_FLOAT*)output;
    if (handle->cfft.len > 0) {
      LC3_cfft_apply(&handle->cfft, dummy, dummy + 1, 2);
    }
    else{
      LC3_iisfft_apply(&handle->iisfft, dummy);
    }
            
    return IIS_FFT_NO_ERROR;
}

IIS_FFT_ERROR LC3_IIS_FFT_Apply_RFFT(HANDLE_IIS_FFT handle, const LC3_FLOAT* in, LC3_FLOAT* out)
{
    if (!handle)
        return IIS_FFT_INTERNAL_ERROR;

    memmove(out, in, sizeof(LC3_FLOAT) * handle->len);

    if (handle->sign == IIS_FFT_BWD)
        LC3_rfft_pre(handle->sine_table, out, handle->len);

    if (handle->cfft.len > 0) {
        LC3_cfft_apply(&handle->cfft, out, out + 1, 2);
    } else {
        LC3_iisfft_apply(&handle->iisfft, out);
    }

    if (handle->sign == IIS_FFT_FWD)
        LC3_rfft_post(handle->sine_table, out, handle->len);

    return IIS_FFT_NO_ERROR;
}
