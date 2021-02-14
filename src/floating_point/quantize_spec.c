/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static LC3_INT sign(LC3_FLOAT x);

LC3_INT sign(LC3_FLOAT x)
{
    if (x > 0)
        return 1;

    if (x < 0)
        return -1;

    return 0;
}

void processQuantizeSpec_fl(LC3_FLOAT x[], LC3_FLOAT gain, LC3_INT xq[], LC3_INT nt, LC3_INT totalBits, LC3_INT* nbits, LC3_INT* nbits2, LC3_INT fs,
                            LC3_INT* lastnzout, LC3_INT* codingdata, LC3_INT* lsbMode, LC3_INT mode, LC3_INT target
#ifdef ENABLE_HR_MODE
                            , LC3_INT hrmode
#endif
)
{

    LC3_INT rateFlag = 0, i = 0, lastnz2 = 0, m = 0, maxlev = 0, k = 0;
    LC3_INT nbits_lsb = 0;
    LC3_INT c         = 0;
    LC3_INT a = 0, b = 0, lev1 = 0, sym = 0, t = 0, pki = 0;
    LC3_INT a1_msb = 0, b1_msb = 0;
    LC3_INT lastnz = 1;
    LC3_FLOAT offset = 0.375;

#ifdef ENABLE_HR_MODE
    if (hrmode)
    {
    	offset = 0.5;
    }
#endif


    /* Quantization */
    for (i = 0; i < nt; i++) {
        xq[i] = trunc(x[i] / gain + offset * sign(x[i]));
#ifdef ENABLE_HR_MODE
        if (hrmode == 0) {
            assert(xq[i] <= 32767 && xq[i] >= -32768);
        }
#else
        assert(xq[i] <= 32767 && xq[i] >= -32768);
#endif
    }

    /* Rate flag */

    if ((fs < 48000 && totalBits > 320 + (fs / 8000 - 2) * 160) || (fs == 48000 && totalBits > 800)) {
        rateFlag = 512;
    }

    /* Init */
    if (mode == 0 && ((fs < 48000 && totalBits >= 640 + (fs / 8000 - 2) * 160) || (fs == 48000 && totalBits >= 1120))) {
        mode = 1;
    }

    /* Last non-uero 2-tuple */

    for (i = nt - 2; i >= 2; i = i - 2) {
        if (xq[i + 1] != 0 || xq[i] != 0) {
            lastnz = i + 1;
            break;
        }
    }


    if (mode < 0) {
        lastnz2 = lastnz + 1;
    } else {
        lastnz2 = 2;
    }

    *nbits  = 0;
    *nbits2 = 0;

    /* Calculate number of estimated bits */

    for (k = 0; k < lastnz; k = k + 2) {
        t = c + rateFlag;
        if (k > nt / 2) {
            t += 256;
        }

        codingdata[0] = t;

        a = abs(xq[k]);
        b = abs(xq[k + 1]);
        m = MAX(a, b);

        if (m == 0) {
            maxlev = -1;
        } else {
            maxlev = floor(LC3_LOG2(MAX(m, 3))) - 1;
        }

        codingdata[1] = maxlev;

        if (mode <= 0) {
            *nbits = *nbits + MIN(a, 1) * 2048;
            *nbits = *nbits + MIN(b, 1) * 2048;
        }

        lev1 = 0;

        while (MAX(a, b) >= 4) {
            pki    = ari_spec_lookup_fl[t + lev1 * 1024];
            *nbits = *nbits + ari_spec_bits_fl[pki][16];

            if (lev1 == 0 && mode > 0) {
                nbits_lsb += 2;
            } else {
                *nbits = *nbits + 2 * 2048;
            }

            a    = a >> 1;
            b    = b >> 1;
            lev1 = MIN(lev1 + 1, 3);
        }

        pki           = ari_spec_lookup_fl[t + lev1 * 1024];
        sym           = a + 4 * b;
        codingdata[2] = sym;
        codingdata += 3;
        *nbits = *nbits + ari_spec_bits_fl[pki][sym];

        if (mode > 0) {
            a1_msb = abs(xq[k]);
            b1_msb = abs(xq[k + 1]);

            if (lev1 > 0) {
                a1_msb = a1_msb >> 1;
                b1_msb = b1_msb >> 1;

                if (a1_msb == 0 && xq[k] != 0) {
                    nbits_lsb++;
                }

                if (b1_msb == 0 && xq[k + 1] != 0) {
                    nbits_lsb++;
                }
            }

            *nbits = *nbits + MIN(a1_msb, 1) * 2048;
            *nbits = *nbits + MIN(b1_msb, 1) * 2048;
        }

        if (mode >= 0 && (abs(xq[k]) != 0 || abs(xq[k + 1]) != 0) && *nbits <= target * 2048) {
            lastnz2 = k + 2;
            *nbits2 = *nbits;
        }

        lev1 = lev1 - 1;
        if (lev1 <= 0) {
            t = 1 + (a + b) * (lev1 + 2);
        } else {
            t = 13 + lev1;
        }

        c = (c & 15) * 16 + t;
    }

    /* Number of bits */
    *nbits = ceil((LC3_FLOAT)*nbits / 2048.0);

    if (mode >= 0) {
        *nbits2 = ceil((LC3_FLOAT)*nbits2 / 2048.0);
    } else {
        *nbits2 = *nbits;
    }

    if (mode > 0) {
        *nbits += nbits_lsb;
        *nbits2 += nbits_lsb;
    }

    /* Truncation of high-frequency coefficients */
    for (i = lastnz2; i <= lastnz; i++) {
        xq[i] = 0;
    }

    /* Truncation of LSBs */
    if (mode > 0 && *nbits > target) {
        *lsbMode = 1;
    } else {
        *lsbMode = 0;
    }

    *lastnzout = lastnz2;
}
