/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processMdctShaping_fl(LC3_FLOAT x[], LC3_FLOAT scf[], const LC3_INT bands_offset[], LC3_INT fdns_npts)
{
    LC3_INT i = 0, j = 0;

    for (i = 0; i < fdns_npts; i++) {
        for (; j < bands_offset[i + 1]; j++) {
            x[j] = x[j] * scf[i];
        }
    }
}
