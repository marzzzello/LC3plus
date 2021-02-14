/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"
#include "fft/iis_fft.c"
#include "fft/iisfft.c"
#include "fft/cfft.c"

void fft_init(Fft* fft, int length)
{
    assert(length % 2 == 0);
    HANDLE_IIS_FFT handle = NULL;
    IIS_FFT_ERROR error = 0;
    
    fft->length = length;
    
    error = LC3_IIS_CFFT_Create(&handle, length, IIS_FFT_FWD);
    
    assert(error == IIS_FFT_NO_ERROR);
    fft->handle = handle;
}

void fft_free(Fft* fft)
{
    IIS_FFT_ERROR error = 0;
    
    if (fft) {
        error = LC3_IIS_CFFT_Destroy((HANDLE_IIS_FFT *) &fft->handle);
        
        assert(error == IIS_FFT_NO_ERROR);
        memset(fft, 0, sizeof(*fft));
    }
}

void fft_apply(Fft* fft, const Complex* input, Complex* output)
{
    IIS_FFT_ERROR error = 0;
    error = LC3_IIS_FFT_Apply_CFFT(fft->handle, input, output);
    
    assert(error == IIS_FFT_NO_ERROR);
}
