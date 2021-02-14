/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processSnsComputeScf_fl(LC3_FLOAT* x, LC3_INT tilt, LC3_INT xLen, LC3_FLOAT* gains, LC3_INT smooth, LC3_FLOAT sns_damping)
{
    LC3_INT   bands_number = 0, d = 0, i = 0, j = 0, n = 0, n2 = 0, n4 = 0, mapping[64] = {0};
    LC3_FLOAT tmp[64] = {0}, x_tmp1[MAX_LEN] = {0}, x_tmp2[MAX_LEN] = {0}, sum = 0, mean = 0, xl4[16] = {0}, nf = 0, xl[64] = {0}, gains_smooth[M] = {0}, ratio = 0;
    LC3_FLOAT W[6] = {1.0 / 12.0, 2.0 / 12.0, 3.0 / 12.0, 3.0 / 12.0, 2.0 / 12.0, 1.0 / 12.0};

    bands_number = xLen;
    assert(bands_number <= 64);

    /* 5 ms */
    if (bands_number < 64) {
        d = 64 - bands_number;

        if (d < xLen)
        {
            j = 0;
            for (i = 0; i < 2 * d; i = i + 2) {
                tmp[i]     = x[j];
                tmp[i + 1] = x[j];
                j++;
            }

            move_float(&tmp[2 * d], &x[d], 64 - 2 * d);
        } else if (ceil(64.0 / (LC3_FLOAT) xLen) == 4)
        {
            ratio = LC3_FABS((LC3_FLOAT) (1.0 - 32.0 / (LC3_FLOAT) xLen));
            n4 = round(ratio * xLen);
            n2 = xLen - n4;
            
            j = 0;
            for(i = 1; i <= n4; i++)
            {
                mapping[j] = i;
                mapping[j + 1] = i;
                mapping[j + 2] = i;
                mapping[j + 3] = i;
                j += 4;
            }
            
            for (i = n4 + 1; i <= n4 + n2; i++)
            {
                mapping[j] = i;
                mapping[j + 1] = i;
                j += 2;
            }
            
            
            for (i = 0; i < 64; i++)
            {
                tmp[i] = x[mapping[i] - 1];
            }
        } else {
            assert(0 && "Unsupported number of bands!");
        }

        move_float(x, tmp, 64);

        bands_number = 64;
        xLen         = bands_number;
    }


    /* Smoothing */

    x_tmp1[0] = x[0];

    move_float(&x_tmp1[1], &x[0], xLen - 1);

    move_float(&x_tmp2[0], &x[1], xLen - 1);

    x_tmp2[xLen - 1] = x[xLen - 1];

    for (i = 0; i < xLen; i++) {
        x[i] = 0.5 * x[i] + 0.25 * x_tmp1[i] + 0.25 * x_tmp2[i];
    }

    /* Pre-emphasis */
    for (i = 0; i < xLen; i++) {
        x[i] = x[i] * LC3_POW(10.0, (LC3_FLOAT)i * (LC3_FLOAT)tilt / ((LC3_FLOAT)bands_number - 1.0) / 10.0);
    }

    /* Noise floor at -40dB */
    for (i = 0; i < 64; i++) {
        sum += x[i];
    }

    mean = sum / (LC3_FLOAT)xLen;

    nf = mean * LC3_POW(10.0, -40.0 / 10.0);
    nf = MAX(nf, LC3_POW(2.0, -32.0));

    for (i = 0; i < 64; i++) {
        if (x[i] < nf) {
            x[i] = nf;
        }
    }

    /* Log-domain */
    for (i = 0; i < 64; i++) {
        xl[i] = LC3_LOG2(x[i]) / 2.0;
    }

    /* Downsampling */
    for (n = 0; n < bands_number / 4; n++) {
        if (n == 0) {
            tmp[0] = xl[0];

            move_float(&tmp[1], &xl[0], 5);

        } else if (n == (bands_number / 4 - 1)) {
            move_float(tmp, &xl[59], 5);

            tmp[5] = xl[63];

        } else {
            move_float(tmp, &xl[n * 4 - 1], ((n * 4 + 5 - 1) - (n * 4 - 1) + 1));
        }

        sum = 0;
        for (i = 0; i < 6; i++) {
            sum += tmp[i] * W[i];
        }

        xl4[n] = sum;
    }


    /* Remove mean and scaling */

    sum = 0;
    for (i = 0; i < bands_number / 4; i++) {
        sum += xl4[i];
    }

    mean = sum / ((LC3_FLOAT)bands_number / 4.0);

    for (i = 0; i < bands_number / 4; i++) {
        gains[i] = sns_damping * (xl4[i] - mean);
    }

    /* Smoothing */
    if (smooth) {
        gains_smooth[0] = (gains[0] + gains[1] + gains[2]) / 3.0;
        gains_smooth[1] = (gains[0] + gains[1] + gains[2] + gains[3]) / 4.0;

        for (i = 2; i < 14; i++) {
            gains_smooth[i] = (gains[i - 2] + gains[i - 1] + gains[i] + gains[i + 1] + gains[i + 2]) / 5.0;
        }

        gains_smooth[M - 2] = (gains[M - 4] + gains[M - 3] + gains[M - 2] + gains[M - 1]) / 4.0;
        gains_smooth[M - 1] = (gains[M - 3] + gains[M - 2] + gains[M - 1]) / 3.0;

        sum = 0;
        for (i = 0; i < M; i++) {
            sum += gains_smooth[i];
        }

        mean = sum / (LC3_FLOAT)M;

        for (i = 0; i < M; i++) {
            gains[i] = 0.5 * (gains_smooth[i] - mean);
        }
    }
}
