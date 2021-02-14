/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static LC3_FLOAT array_max_abs(LC3_FLOAT* x, LC3_INT len);

LC3_FLOAT array_max_abs(LC3_FLOAT* x, LC3_INT len)
{
    LC3_FLOAT max = 0, tmp_abs = 0;
    LC3_INT   i = 0;

    for (i = 0; i < len; i++) {
        tmp_abs = LC3_FABS(x[i]);
        if (tmp_abs > max) {
            max = tmp_abs;
        }
    }

    return max;
}

void processEstimateGlobalGain_fl(LC3_FLOAT x[], LC3_INT lg, LC3_INT nbitsSQ, LC3_FLOAT* gain, LC3_INT* quantizedGain,
                                  LC3_INT* quantizedGainMin, LC3_INT quantizedGainOff, LC3_FLOAT* targetBitsOff,
                                  LC3_INT* old_targetBits, LC3_INT old_specBits
#ifdef ENABLE_HR_MODE
                                  , LC3_INT hrmode , LC3_INT regBits, LC3_FLOAT frame_ms
#endif
)
{

    LC3_INT   i = 0, N = 0, offset = 0, j = 0, iszero = 0;
    LC3_FLOAT g_min = 0, x_max = 0, tmp = 0, ind = 0, ind_min = 0, target = 0, fac = 0, ener = 0;
    LC3_FLOAT en[MAX_LEN / 4] = {0};
    LC3_FLOAT reg_val = 0;

    if (*old_targetBits < 0) {
        *targetBitsOff = 0;
    } else {
        tmp            = MIN(40, MAX(-40, *targetBitsOff + *old_targetBits - old_specBits));
        *targetBitsOff = 0.8 * *targetBitsOff + 0.2 * tmp;
    }

    *old_targetBits = nbitsSQ;
    nbitsSQ         = nbitsSQ + round(*targetBitsOff);

    x_max = array_max_abs(x, lg);

#ifdef ENABLE_HR_MODE
    if (hrmode && regBits > 0)
    {
        LC3_FLOAT M0 = 1e-5, M1 = 1e-5, rB_offset;
        LC3_FLOAT thresh = 2*frame_ms;
        for (i = 0; i < lg; i++)
        {
        		M0 += fabs(x[i]);
        		M1 += i*fabs(x[i]);
        }

        rB_offset = 8 * (1 - MIN(M1/M0, thresh) / thresh);
        reg_val = x_max * LC3_POW(2,-regBits - rB_offset);
    }
#endif

    if (x_max == 0) {
        ind_min         = quantizedGainOff;
        ind             = 0;
        *old_targetBits = -1;
    } else {
#ifdef ENABLE_HR_MODE
        if (hrmode == 1) {
        	g_min = x_max / (32768 * 256 - 2);
        } else {
            g_min = x_max / (32768 - 0.375);
        }
#else
        g_min = x_max / (32768 - 0.375);
#endif
        ind_min = ceil(28.0 * LC3_LOG10(g_min));
        assert(ind_min <= (255 + quantizedGainOff));

        N = lg;

        j = 0;
        for (i = 0; i < N; i = i + 4) {
            tmp = x[i] * x[i];
            tmp += x[i + 1] * x[i + 1];
            tmp += x[i + 2] * x[i + 2];
            tmp += x[i + 3] * x[i + 3];
            en[j] = (28.0 / 20.0) * (7 + 10.0 * LC3_LOG10(tmp + reg_val + LC3_POW(2, -31)));
            j++;
        }

        target = (28.0 / 20.0) * (1.4) * nbitsSQ;
        fac    = 256;
        offset = 255 + quantizedGainOff;

        for (i = 0; i < 8; i++) {
            fac    = fac * 0.5;
            offset = offset - fac;
            ener   = 0;
            iszero = 1;

            for (j = N / 4 - 1; j >= 0; j--) {
                tmp = en[j] - offset;

                if (tmp < (7.0) * (28.0 / 20.0)) {
                    if (iszero == 0) {
                        ener = ener + (2.7) * (28.0 / 20.0);
                    }
                } else {
                    if (tmp > (50.0) * (28.0 / 20.0)) {
                        ener = ener + 2.0 * tmp - (50.0) * (28.0 / 20.0);
                    } else {
                        ener = ener + tmp;
                    }

                    iszero = 0;
                }
            }

            if (ener > target && iszero == 0) {
                offset = offset + fac;
            }
        }

        if (offset < ind_min) {
            *old_targetBits = -1;
        }

        ind = MAX(ind_min, offset) - quantizedGainOff;
    }

    *quantizedGainMin = ind_min;
    *quantizedGain    = ind;

    *gain = LC3_POW(10.0, ((ind + quantizedGainOff) / 28.0));
}
