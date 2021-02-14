/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "../functions.h"

#ifndef CFFT_H
#define CFFT_H

typedef struct
{
    LC3_INT len;
    LC3_INT sign;
    LC3_FLOAT* table;
} Cfft;

/* macro to check if cfft supports len  */
#define CFFT_IS_POWER_OF_TWO(n) ((n != 0) && ((n & (~n + 1)) == n))
#define CFFT_SUPPORT(len) (len >= 4 && len <= 1024 && CFFT_IS_POWER_OF_TWO(len))
#define CFFT_PLAN_SUPPORT(len) (len >= 4 && CFFT_IS_POWER_OF_TWO(len))

/**
 * \brief fft_radix2
          The function serves as wrapper for the forward and inverse complex fft.

 * \param[i/o] re:          real input / real output
 * \param[i/o] im:          imag input / imag output
 * \param[i  ] sizeOfFft:   size of fft
 * \param[i  ] s:           stride of real and imag input / output
 * \param[i  ] iSign:       forward fft: -1 / inverse fft: 1

 * \return none
 */

void LC3_cfft(LC3_FLOAT* re, LC3_FLOAT* im, LC3_INT sizeOfFft, LC3_INT stride, LC3_INT sign);

LC3_INT LC3_cfft_plan(Cfft* handle, LC3_INT length, LC3_INT sign);
void LC3_cfft_apply(Cfft* handle, LC3_FLOAT* re, LC3_FLOAT* im, LC3_INT stride);
void LC3_cfft_free(Cfft* handle);

#endif /* FFT_RADIX2_H */
