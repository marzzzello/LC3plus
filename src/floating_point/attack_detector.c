/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void attack_detector_fl(LC3_FLOAT* in, LC3_INT frame_size, LC3_INT* lastAttackPosition, LC3_FLOAT* accNrg, LC3_INT* attackFlag,
                        LC3_FLOAT* attdec_filter_mem, LC3_INT attackHandlingOn)
{
    LC3_FLOAT f_sig[160] = {0}, block_nrg[4] = {0}, sum = 0, tmpEne = 0, *ptr = NULL, tmp[162] = {0};
    LC3_INT   i = 0, j = 0, attackPosition = 0;
    LC3_FLOAT mval = 0;

    ptr = &tmp[2];

    if (attackHandlingOn) {
        /* Decimate 96, 48 and 32 kHz signals to 16 kHz */
        if (frame_size == 960) {
            for (i = 0; i < frame_size;) {
                ptr[j] = in[i] + in[i + 1] + in[i + 2] + in[i + 3] + in[i + 4] + in[i + 5];
                i      = i + 6;
                j++;
            }
            mval = 1e-5;
        } else if (frame_size == 480) {
            j = 0;
            for (i = 0; i < frame_size;) {
                ptr[j] = (in[i] + in[i + 1] + in[i + 2]);
                i      = i + 3;
                j++;
            }
        } else if (frame_size == 320) {
            j = 0;
            for (i = 0; i < frame_size;) {
                ptr[j] = (in[i] + in[i + 1]);
                i      = i + 2;
                j++;
            }
        } else if (frame_size == 240) {
            j = 0;
            for (i = 0; i < frame_size;) {
                ptr[j] = (in[i] + (in[i + 1] + in[i + 2]) / 2.0);
                i      = i + 3;
                j++;
            }
        }

        /* Filter */
        ptr[-2] = (LC3_FLOAT)attdec_filter_mem[0];
        ptr[-1] = (LC3_FLOAT)attdec_filter_mem[1];

        attdec_filter_mem[0] = ptr[158];
        attdec_filter_mem[1] = ptr[159];

        for (i = 159; i >= 0; i--) {
            tmpEne = 0;

            tmpEne += ptr[i] * 0.375;
            tmpEne += ptr[i - 1] * (-0.5);
            tmpEne += ptr[i - 2] * (0.125);

            f_sig[i] = tmpEne;
        }

        for (i = 0; i < 4; i++) {
            sum = 0;
            for (j = 0; j < 40; j++) {
                sum += f_sig[j + i * 40] * f_sig[j + i * 40];
            }

            block_nrg[i] = sum;
        }

        *attackFlag    = 0;
        attackPosition = -1;

        for (i = 0; i < 4; i++) {
            tmpEne = block_nrg[i] / 8.5;

            if (tmpEne > MAX(*accNrg, mval)) {
                *attackFlag    = 1;
                attackPosition = i + 1;
            }

            *accNrg = MAX(block_nrg[i], 0.25 * (*accNrg));
        }

        if (*lastAttackPosition > 2) {
            *attackFlag = 1;
        }

        *lastAttackPosition = attackPosition;
    }
}
