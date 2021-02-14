/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static void filter_olpa(LC3_FLOAT* in, LC3_FLOAT* out, const LC3_FLOAT* buf, LC3_FLOAT len_buf, LC3_INT len_input);
static LC3_INT  searchMaxIndice(LC3_FLOAT* in, LC3_INT len);

void filter_olpa(LC3_FLOAT* in, LC3_FLOAT* out, const LC3_FLOAT* buf, LC3_FLOAT len_buf, LC3_INT len_input)
{
    LC3_INT   i = 0, j = 0;
    LC3_FLOAT sum = 0;
    /* a = 1, so denominator == 1, nothing to do here */

    for (i = 0; i < len_input; i++) {
        j   = 0;
        sum = 0;
        for (j = 0; (j < len_buf) && (j <= i); j++) {
            sum += buf[j] * in[i - j];
        }

        out[i] = sum;
    }
}

LC3_INT searchMaxIndice(LC3_FLOAT* in, LC3_INT len)
{
    LC3_INT   max_i = 0, i = 0;
    LC3_FLOAT max = in[0];

    if (len <= 0) {
        return -128;
    }

    for (i = 0; i < len; i++) {
        if (in[i] > max) {
            max   = in[i];
            max_i = i;
        }
    }

    return max_i;
}

void processOlpa_fl(LC3_FLOAT* wsp, LC3_FLOAT* mem_lp_decim2, LC3_FLOAT* mem_old_d_wsp, LC3_INT* mem_old_T0, LC3_INT* T0_out,
                    LC3_FLOAT* normcorr_out, LC3_INT wsp_len)
{
    LC3_FLOAT norm_corr = 0, sum = 0, sum_sq1 = 0, sum_sq2 = 0, tmp[480] = {0}, sum_tmp = 0, buf_tmp[480] = {0},
          norm_corr2 = 0;
    LC3_FLOAT buf[178] = {0}, filt_out[131] = {0}, d_wsp[64] = {0}, R0[98] = {0}, R[98] = {0}; /* constant length */
    LC3_INT   i = 0, j = 0, N = 0, d_N = 0, buf_len = 0, L_min = 0, L_max = 0, T0 = 0, T02 = 0, L_min2 = 0, L_max2 = 0, L = 0;

    N       = wsp_len;
    d_N     = N / 2;
    buf_len = wsp_len + 3;

    /* Downsampling */

    move_float(buf, mem_lp_decim2, 3);

    move_float(&buf[3], wsp, wsp_len);

    move_float(mem_lp_decim2, &buf[N], buf_len - N);

    filter_olpa(buf, filt_out, olpa_down2, 5, buf_len);

    j = 0;
    for (i = 4; i < buf_len; i = i + 2) {
        d_wsp[j] = filt_out[i];
        j++;
    }

    /* Correlation */

    move_float(buf, mem_old_d_wsp, 114);

    move_float(&buf[114], d_wsp, 64);

    move_float(mem_old_d_wsp, &buf[d_N], 178 - d_N);

    L_min = 17;
    L_max = 114;

    for (i = L_min; i <= L_max; i++) {
        sum = 0;

        move_float(buf_tmp, &buf[L_max - i], d_N);

        for (j = 0; j < d_N; j++) {
            sum += d_wsp[j] * buf_tmp[j];
        }

        R0[i - L_min] = sum;
    }

    move_float(R, R0, 98);

    for (i = 0; i < 98; i++) {
        R0[i] = R0[i] * olpa_acw[i];
    }

    L  = searchMaxIndice(R0, 98);
    T0 = L + L_min;

    move_float(tmp, &buf[L_max - T0], d_N);

    sum_tmp = 0;
    sum_sq1 = 0;
    sum_sq2 = 0;

    for (i = 0; i < d_N; i++) {
        sum_tmp += d_wsp[i] * tmp[i];
        sum_sq1 += d_wsp[i] * d_wsp[i];
        sum_sq2 += tmp[i] * tmp[i];
    }

    sum_sq1 = sum_sq1 * sum_sq2;
    sum_sq1 = LC3_SQRT(sum_sq1) + LC3_POW(10.0, -5.0);

    norm_corr = sum_tmp / sum_sq1;
    norm_corr = MAX(0, norm_corr);

    L_min2 = MAX(L_min, *mem_old_T0 - 4);
    L_max2 = MIN(L_max, *mem_old_T0 + 4);

    /* Distance can be negative here */
    if (((L_max2 - L_min + 1) - (L_min2 - L_min)) > 0) {
        move_float(buf_tmp, &R[L_min2 - L_min], ((L_max2 - L_min + 1) - (L_min2 - L_min)));
    }

    L = searchMaxIndice(buf_tmp, (L_max2 - L_min + 1) - (L_min2 - L_min));

    if (L == -128) {
        T02 = -128;
    } else {
        T02 = L + L_min2;
    }

    if ((T02 != T0) && (T02 != -128)) {
        move_float(tmp, &buf[L_max - T02], d_N);

        sum_tmp = 0;
        sum_sq1 = 0;
        sum_sq2 = 0;
        for (i = 0; i < d_N; i++) {
            sum_tmp += d_wsp[i] * tmp[i];
            sum_sq1 += d_wsp[i] * d_wsp[i];
            sum_sq2 += tmp[i] * tmp[i];
        }

        sum_sq1 = sum_sq1 * sum_sq2;
        sum_sq1 = LC3_SQRT(sum_sq1) + LC3_POW(10.0, -5.0);

        norm_corr2 = sum_tmp / sum_sq1;

        norm_corr2 = MAX(0, norm_corr2);

        if (norm_corr2 > (norm_corr * 0.85)) {
            T0        = T02;
            norm_corr = norm_corr2;
        }
    }

    *mem_old_T0 = T0;

    T0 = T0 * 2.0;

    *T0_out       = T0;
    *normcorr_out = norm_corr;
}
