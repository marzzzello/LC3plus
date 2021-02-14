/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void process_cutoff_bandwidth(LC3_FLOAT *d_fl, LC3_INT len, LC3_INT bw_bin)
{
    LC3_INT i = 0;
    
    if (len > bw_bin){
        for (i = -1; i < 3; i++) {
            d_fl[bw_bin + i] = d_fl[bw_bin + i] * LC3_POW(2, -(i + 2));
        }

        for (i = bw_bin + 3; i < len; i++) {
            d_fl[i] = 0;
        } 
    }
}
