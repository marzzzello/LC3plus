/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processSnsInterpolateScf_fl(LC3_FLOAT* gains, LC3_INT encoder_side, LC3_INT bands_number, LC3_FLOAT* gains_int)
{
    LC3_INT   i = 0, n = 0, d = 0, n4 = 0;
    LC3_FLOAT tmp[MAX_BANDS_NUMBER_PLC] = {0}, ratio = 0;

    /* Interpolation */

    gains_int[0] = gains[0];
    gains_int[1] = gains[0];

    for (n = 0; n <= 14; n++) {
        gains_int[n * 4 + 2] = gains[n] + (gains[n + 1] - gains[n]) / 8.0;
        gains_int[n * 4 + 3] = gains[n] + 3.0 * (gains[n + 1] - gains[n]) / 8.0;
        gains_int[n * 4 + 4] = gains[n] + 5.0 * (gains[n + 1] - gains[n]) / 8.0;
        gains_int[n * 4 + 5] = gains[n] + 7.0 * (gains[n + 1] - gains[n]) / 8.0;
    }

    gains_int[62] = gains[15] + (gains[15] - gains[14]) / 8.0;
    gains_int[63] = gains[15] + 3.0 * (gains[15] - gains[14]) / 8.0;

    /* For 5ms */

    if (bands_number < 64) {
        d = 64 - bands_number;

        if (d < 32)
        {
            i = 0;
            for (n = 0; n < 2 * d; n = n + 2) {
                tmp[i] = (gains_int[n] + gains_int[n + 1]) / 2.0;
                i++;
            }

            for (n = 1; n < d; n++) {
                gains_int[n] = gains_int[2 * n];
            }

            for (n = 2 * d; n < 64; n++) {
                gains_int[n - d] = gains_int[n];
            }

            move_float(gains_int, tmp, d);
        } else if (ceil(64.0 / (LC3_FLOAT) bands_number) == 4)
        {
            ratio = LC3_FABS((LC3_FLOAT) (1.0 - 32.0 / (LC3_FLOAT) bands_number));
            n4 = round(ratio * bands_number);
            
            for (i = 0; i < n4; i++)
            {
                tmp[i] = (gains_int[4 * i] + gains_int[4 * i + 1] + gains_int[4 * i + 2] + gains_int[4 * i + 3]) / 4.0;
            }
            
            for (i = 0; i < bands_number - n4; i++)
            {
                tmp[n4 + i] = (gains_int[4 * n4 + 2 * i] + gains_int[4 * n4 + 2 * i + 1]) / 2.0;
            }
            
            move_float(gains_int, tmp, bands_number);
        } else {
            assert(0 && "Unsupported number of bands!");
        }
    }

    /* Inversion at encoder-side */

    if (encoder_side == 1) {
        for (n = 0; n < bands_number; n++) {
            gains_int[n] = -gains_int[n];
        }
    }

    /* Linear domain */

    for (n = 0; n < bands_number; n++) {
        gains_int[n] = LC3_POW(2, gains_int[n]);
    }
}
