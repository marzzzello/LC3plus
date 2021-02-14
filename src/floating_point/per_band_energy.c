/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processPerBandEnergy_fl(LC3_INT bands_number, const LC3_INT* acc_coeff_per_band, LC3_FLOAT* d2, LC3_FLOAT* d)
{
    LC3_INT   i = 0, j = 0, start = 0, stop = 0;
    LC3_FLOAT sum = 0;

    for (i = 0; i < bands_number; i++) {
        sum   = 0;
        start = acc_coeff_per_band[i];
        stop  = acc_coeff_per_band[i + 1];

        for (j = start; j < stop; j++) {
            sum += d[j] * d[j];
        }

        sum   = sum / (LC3_FLOAT)(stop - start);
        d2[i] = sum;
    }
}
