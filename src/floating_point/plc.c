/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processPlcApply(LC3_FLOAT* spec_prev, LC3_INT L_spec, LC3_INT* nbLostCmpt, LC3_FLOAT* cum_alpha, LC3_INT* seed, LC3_FLOAT* spec_out)
{
    LC3_INT   i     = 0;
    LC3_FLOAT alpha = 0;

    *nbLostCmpt = *nbLostCmpt + 1;

    /* Get damping factor */
    if (*nbLostCmpt < 4) {
        alpha = 1.0;
    } else if (*nbLostCmpt < 8) {
        alpha = 0.9;
    } else {
        alpha = 0.85;
    }

    *cum_alpha = *cum_alpha * alpha;

    /* Noise Substitution */
    for (i = 0; i < L_spec; i++) {
        *seed = 16831 + *seed * 12821;
        *seed = *seed - round(*seed * LC3_POW(2, -16)) * LC3_POW(2, 16);

        if (*seed == 32768) {
            *seed = *seed - 32768;
        }

        if (*seed < 0) {
            spec_out[i] = spec_prev[i] * (-(*cum_alpha));
        } else {
            spec_out[i] = spec_prev[i] * (*cum_alpha);
        }
    }
}

void processPlcUpdate(LC3_FLOAT* spec_cur, LC3_FLOAT* spec_prev, LC3_INT len, LC3_INT* nbLostCmpt, LC3_FLOAT* cum_alpha)
{
    move_float(spec_prev, spec_cur, len);
    *nbLostCmpt = 0;
    *cum_alpha  = 1;
}

