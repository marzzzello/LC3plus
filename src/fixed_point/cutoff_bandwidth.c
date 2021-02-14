/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"


void process_cutoff_bandwidth(Word32 d_fx[], Word16 len, Word16 bw_bin)
{
    Counter i = 0;
    if (len > bw_bin){
        /* roll off */
        for (i = -1; i < 3; i++) {
            d_fx[bw_bin + i] = L_shr(d_fx[bw_bin + i], add(i, 2));
        }

        for (i = bw_bin + 3; i < len; i++) {
            d_fx[i] = 0; move32();
        } 
    }
}

