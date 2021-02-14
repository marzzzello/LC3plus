/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void dct2_init(Dct2* dct, int length)
{
    assert(length % 2 == 0);
    assert(length <= MAX_LEN);
    dct->length = length;
    fft_init(&dct->fft, length);
}

void dct2_free(Dct2* dct)
{
    if (dct) {
        fft_free(&dct->fft);
        memset(dct, 0, sizeof(*dct));
    }
}

void dct2_apply(Dct2* dct, const LC3_FLOAT* input, LC3_FLOAT* output)
{
    assert(input != output);
    Complex   tmp1[MAX_LEN];
    Complex   tmp2[MAX_LEN];
    int       i   = 0;
    const int len = dct->length;

    for (i = 0; i < len / 2; i++) {
        tmp1[i]           = cmplx(input[i * 2], 0);
        tmp1[len - i - 1] = cmplx(input[i * 2 + 1], 0);
    }

    fft_apply(&dct->fft, tmp1, tmp2);

    for (i = 0; i < len; i++) {
        Complex c = cmul(cexpi(-M_PI * i / (2 * len)), cmplx(2 / LC3_SQRT(2 * len), 0));
        output[i] = cmul(tmp2[i], c).r;
    }
    output[0] /= LC3_SQRT(2);
}

void dct3_init(Dct3* dct, int length)
{
    assert(length % 2 == 0);
    assert(length <= MAX_LEN);
    dct->length = length;
    fft_init(&dct->fft, length);
}

void dct3_free(Dct3* dct)
{
    if (dct) {
        fft_free(&dct->fft);
        memset(dct, 0, sizeof(*dct));
    }
}

void dct3_apply(Dct3* dct, const LC3_FLOAT* input, LC3_FLOAT* output)
{
    assert(input != output);
    Complex     tmp1[MAX_LEN] = {{0}};
    Complex     tmp2[MAX_LEN] = {{0}};
    int         i    = 0;
    const int   len  = dct->length;
    const LC3_FLOAT norm = 1.0 / len;

    for (i = 0; i < len; i++) {
        Complex c = cmul(cexpi(M_PI * i / (2 * len)), cmplx(LC3_SQRT(2 * len), 0));
        c         = cmul(c, cmplx(input[i], 0));
        tmp1[i]   = cmplx(c.i, c.r);
    }
    tmp1[0].i /= LC3_SQRT(2);

    fft_apply(&dct->fft, tmp1, tmp2);

    for (i = 0; i < len / 2; i++) {
        output[i * 2]     = tmp2[i].i * norm;
        output[i * 2 + 1] = tmp2[len - i - 1].i * norm;
    }
}

void dct4_init(Dct4* dct, int length)
{
    assert(length % 2 == 0);
    assert(length <= MAX_LEN);
    int i       = 0;
    dct->length = length;
    dct->twid1  = calloc(sizeof(*dct->twid1), length / 2);
    dct->twid2  = calloc(sizeof(*dct->twid2), length / 2);
    for (i = 0; i < length / 2; i++) {
        dct->twid1[i] = cexpi(-M_PI * (i + 0.25) / length);
        dct->twid2[i] = cexpi(-M_PI * i / length);
    }
    fft_init(&dct->fft, length / 2);
}

void dct4_free(Dct4* dct)
{
    if (dct) {
        free(dct->twid1);
        free(dct->twid2);
        fft_free(&dct->fft);
        memset(dct, 0, sizeof(*dct));
    }
}

void dct4_apply(Dct4* dct, const LC3_FLOAT* input, LC3_FLOAT* output)
{
    assert(input != output);
    Complex     tmp2[MAX_LEN / 2];
    int         i    = 0;
    Complex*    tmp1 = (Complex*)output;
    const int   len  = dct->length;
    const LC3_FLOAT norm = 1.0 / LC3_SQRT(len / 2);

    for (i = 0; i < len / 2; i++) {
        tmp1[i] = cmul(cmplx(input[i * 2], input[len - i * 2 - 1]), dct->twid1[i]);
    }

    fft_apply(&dct->fft, tmp1, tmp2);

    for (i = 0; i < len / 2; i++) {
        Complex t               = cmul(tmp2[i], dct->twid2[i]);
        output[i * 2]           = t.r * norm;
        output[len - i * 2 - 1] = -t.i * norm;
    }
}
