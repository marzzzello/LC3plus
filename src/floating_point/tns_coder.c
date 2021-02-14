/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static void xcorr(LC3_FLOAT* in, LC3_FLOAT* out, LC3_INT lag, LC3_INT inLen);
static void levdown(LC3_FLOAT* anxt, LC3_FLOAT* out_a, LC3_INT* len);
static void poly2rc(LC3_FLOAT* a, LC3_FLOAT* out, LC3_INT len);
static LC3_INT  findRC_idx(const LC3_FLOAT* in1, const LC3_FLOAT* in2, LC3_FLOAT checkValue);

void xcorr(LC3_FLOAT* in, LC3_FLOAT* out, LC3_INT lag, LC3_INT inLen)
{
    LC3_INT   i = 0, m = 0;
    LC3_FLOAT sum = 0, tmp_buf[MAX_LEN] = {0};

    for (m = -lag; m <= lag; m++) {
        /* Append zeros and input vector */

        zero_float(tmp_buf, abs(m));

        move_float(&tmp_buf[abs(m)], in, inLen - abs(m));

        /* Calculate sum */
        sum = 0;

        for (i = 0; i < inLen; i++) {
            sum += in[i] * tmp_buf[i];
        }

        out[m + lag] = sum;
    }
}

void levinsonDurbin(LC3_FLOAT* r, LC3_FLOAT* out_lev, LC3_FLOAT* rc_unq, LC3_FLOAT* error, LC3_INT len)
{
    LC3_INT   t = 0, i = 0, j = 0;
    LC3_FLOAT g = 0, v = 0, sum = 0, buf_tmp[MAX_LEN] = {0};

    g          = r[1] / r[0];
    out_lev[0] = g;

    v         = (1.0 - g * g) * r[0];
    rc_unq[0] = -g;

    for (t = 1; t < len; t++) {
        zero_float(buf_tmp, len + 1);

        sum = 0;
        for (i = 1; i <= t; i++) {
            sum += out_lev[i - 1] * r[i];
        }

        g = (r[t + 1] - sum) / v;

        j = 1;
        for (i = t - 1; i >= 0; i--) {
            buf_tmp[j] = out_lev[j - 1] - g * out_lev[i];
            j++;
        }

        move_float(&out_lev[1], &buf_tmp[1], len);

        out_lev[0] = g;

        v         = v * (1 - g * g);
        rc_unq[t] = -g;
    }

    /* Reorder out_lev */
    out_lev[0] = 1;
    j          = 1;
    for (i = len - 1; i >= 0; i--) {
        buf_tmp[j] = -out_lev[i];
        j++;
    }

    move_float(&out_lev[1], &buf_tmp[1], (len - 1));

    out_lev[len] = rc_unq[len - 1];

    *error = v;
}

void levdown(LC3_FLOAT* anxt, LC3_FLOAT* out_a, LC3_INT* len)
{
    LC3_INT   i = 0, j = 0;
    LC3_FLOAT tmp_buf[8] = {0}, tmp_buf1[8] = {0}, tmp_buf2[8] = {0}, knxt = 0;

    /* Initial length = 9 */

    /* Drop the leading 1 */

    move_float(&tmp_buf[0], &anxt[1], (*len - 1));

    *len = *len - 1; /* Lenght = 8 */

    /* Last coefficient */
    knxt = tmp_buf[*len - 1]; /* At [7] */

    *len = *len - 1; /* Lenght = 7 */

    move_float(tmp_buf1, tmp_buf, *len);

    j = 0;
    for (i = *len - 1; i >= 0; i--) {
        tmp_buf2[j] = knxt * tmp_buf[i];
        j++;
    }

    out_a[0] = 1;
    for (i = 0; i < *len; i++) {
        out_a[i + 1] = (tmp_buf1[i] - tmp_buf2[i]) / (1.0 - (LC3_FABS(knxt)) * (LC3_FABS(knxt)));
    }

    *len = *len + 1; /* Length = 8 */
}

void poly2rc(LC3_FLOAT* a, LC3_FLOAT* out, LC3_INT len)
{
    LC3_INT   k = 0, i = 0, len_old = 0;
    LC3_FLOAT buf[9] = {0};

    len_old = len;

    zero_float(out, len - 1);

    /* Length = 9 */

    /* Normalize */
    for (i = 0; i < len; i++) {
        a[i] = a[i] / a[0];
    }

    out[len - 1] = a[len - 1];

    /* Process */
    for (k = len - 2; k >= 0; k--) {
        levdown(a, buf, &len);
        out[k] = buf[len - 1]; /* Store last value */

        move_float(a, buf, len);
    }

    /* Shift output array by one to the left to lose leading 1 */
    for (i = 0; i < len_old - 1; i++) {
        out[i] = out[i + 1];
    }
}

LC3_INT findRC_idx(const LC3_FLOAT* in1, const LC3_FLOAT* in2, LC3_FLOAT checkValue)
{
    LC3_INT i = 0, ret = 0;

    for (i = 0; i < 17; i++) {
        if (checkValue <= in1[i] && checkValue > in2[i]) {
            ret = i;
        }
    }

    return ret;
}

void processTnsCoder_fl(LC3_FLOAT* x, LC3_INT bw_cutoff_idx, LC3_INT bw_fcbin, LC3_INT fs, LC3_INT N, LC3_INT frame_dms, LC3_INT nBits,
                        LC3_INT* order_out, LC3_INT* rc_idx, LC3_INT* tns_numfilters, LC3_INT* bits_out)
{
    LC3_INT i = 0, stopfreq[2] = {0}, startfreq[2] = {0}, f = 0, numfilters = 0, maxOrder = 0, bits = 0, sub = 0,
        subdiv_startfreq = 0, subdiv_stopfreq = 0, j = 0, rc_idx_tmp[8] = {0}, order_tmp[8] = {0}, tmp = 0, tns = 0;
    LC3_FLOAT minPGfac = 0, minPredictionGain = 0, maxPG = 0, xcorr_out[MAX_LEN] = {0}, buf_tmp[MAX_LEN] = {0}, sum = 0,
          subdiv_len = 0, nSubdivisions = 0, r[9] = {0}, out_lev[9] = {0}, rc_unq[9] = {0}, error_lev = 0, predGain = 0,
          alpha = 0, rc[8] = {0}, st[9] = {0}, s = 0, tmpSave = 0, tmp_fl = 0;
    const LC3_INT* order;

    /* Init */

    if (fs >= 32000 && frame_dms >= 50) {
        numfilters = 2;
    } else {
        numfilters = 1;
    }

    if (N > 40 * ((LC3_FLOAT) (frame_dms) / 10.0)) {
        N  = 40 * ((LC3_FLOAT) (frame_dms) / 10.0);
        fs = 40000;
    }

    if (numfilters == 1) {
        startfreq[0] = floor(600 * N * 2 / fs) + 1;
        stopfreq[0]  = N;
    } else {
        startfreq[0] = floor(600 * N * 2 / fs) + 1;
        startfreq[1] = N / 2 + 1;
        stopfreq[0]  = N / 2;
        stopfreq[1]  = N;
    }
    
    switch (frame_dms)
    {
        case 25:
            maxOrder      = 4;
            nSubdivisions = 2.0;
            break;
        case 50:
            maxOrder      = 4;
            nSubdivisions = 2.0;
            break;
        case 100:
            maxOrder      = 8;
            nSubdivisions = 3.0;
            break;
    }

    minPGfac          = 0.85;
    maxPG             = 2;
    minPredictionGain = 1.5;

    if ((frame_dms >= 50 && nBits >= 48 * ((LC3_FLOAT) frame_dms / 10.0)) || frame_dms == 25) {
        maxPG = minPredictionGain;
    }

    if ((frame_dms >= 50 && nBits >= 48 * ((LC3_FLOAT) frame_dms / 10.0)) || frame_dms == 25) {
        order = order1_tns;
    } else {
        order = order2_tns;
    }
    
    /* Processing */
    if (bw_cutoff_idx >= 3 && numfilters == 2) {
        numfilters   = 2;
        startfreq[1] = bw_fcbin / 2 + 1;
        stopfreq[0]  = bw_fcbin / 2;
        stopfreq[1]  = bw_fcbin;
    } else {
        numfilters  = 1;    
        stopfreq[0] = bw_fcbin;
    }

    bits = 0;

    for (f = 0; f < numfilters; f++) {
        subdiv_len = ((LC3_FLOAT)stopfreq[f] + 1.0 - (LC3_FLOAT)startfreq[f]) / nSubdivisions;

        zero_float(r, 9);

        for (sub = 1; sub <= nSubdivisions; sub++) {
            subdiv_startfreq = floor(subdiv_len * (sub - 1)) + startfreq[f] - 1;
            subdiv_stopfreq  = floor(subdiv_len * sub) + startfreq[f] - 1;

            sum = 0;
            for (i = subdiv_startfreq; i < subdiv_stopfreq; i++) {
                sum += x[i] * x[i];
            }

            if (sum == 0) {
                zero_float(r, 9);
                r[0] = 1;
                break;
            }

            move_float(buf_tmp, &x[subdiv_startfreq], subdiv_stopfreq - subdiv_startfreq);

            xcorr(buf_tmp, xcorr_out, maxOrder, subdiv_stopfreq - subdiv_startfreq);

            j = 0;
            for (i = maxOrder; i >= 0; i--) {
                r[j] = r[j] + xcorr_out[i] / sum;
                j++;
            }
        }

        for (i = 0; i <= maxOrder; i++) {
            r[i] = r[i] * lagw_tns[i];
        }

        levinsonDurbin(r, out_lev, rc_unq, &error_lev, maxOrder);

        predGain = r[0] / error_lev;

        if (predGain > minPredictionGain) {
            tns = 1;
        } else {
            tns = 0;
        }

        bits++;

        if (tns == 1) {
            /* LPC weighting */
            if (predGain < maxPG) {
                alpha = (maxPG - predGain) * (minPGfac - 1.0) / (maxPG - minPredictionGain) + 1.0;

                for (i = 0; i <= maxOrder; i++) {
                    out_lev[i] = out_lev[i] * LC3_POW(alpha, i);
                }

                poly2rc(out_lev, rc_unq, maxOrder + 1);
            }

            /* PARCOR Quantization */
            for (i = 0; i < maxOrder; i++) {
                rc_idx_tmp[i] = findRC_idx(&quants_thr_tns[1], &quants_thr_tns[0], rc_unq[i]);
            }
            
            /* Filter Order */
            j = 0;
            for (i = 0; i < maxOrder; i++) {
                rc[i] = quants_pts_tns[rc_idx_tmp[i]];

                if (rc[i] != 0) {
                    order_tmp[j] = i + 1;
                    j++;
                }
            }

            order_out[f] = order_tmp[j - 1];
            tmp = order[order_out[f] - 1];

            /* Huffman Coding of PARCOR coefficients */
            for (i = 0; i <= order_out[f] - 1; i++) {
                tmp += huff_bits_tns[i][rc_idx_tmp[i]];
            }

            bits = bits + ceil((LC3_FLOAT)tmp / 2048.0);

            j = 0;
            for (i = f * 8; i <= f * 8 + order_out[f] - 1; i++) {
                rc_idx[i] = rc_idx_tmp[j];
                j++;
            }
        }

        /* Filtering */
        if (tns == 1) {
            for (i = startfreq[f]; i <= stopfreq[f]; i++) {
                s       = x[i - 1];
                tmpSave = s;

                for (j = 0; j < order_out[f] - 1; j++) {
                    tmp_fl = rc[j] * s + st[j];
                    s += rc[j] * st[j];

                    st[j]   = tmpSave;
                    tmpSave = tmp_fl;
                }

                s += rc[order_out[f] - 1] * st[order_out[f] - 1];

                st[order_out[f] - 1] = tmpSave;
                x[i - 1]             = s;
            }
        }
    }

    *tns_numfilters = numfilters;
    *bits_out       = bits;
}
