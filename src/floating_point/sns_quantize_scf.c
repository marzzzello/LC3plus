/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static void pvq_dec(LC3_INT k, LC3_INT m, LC3_INT LS_ind, LC3_INT MPVQ_ind, LC3_INT* pulses);
static LC3_INT  find_last_indice_le(LC3_INT compare, const LC3_INT* array, LC3_INT len);
static void pvq_enc(LC3_INT* pulses, LC3_INT* LS_ind, LC3_INT* MPVQ_ind, LC3_INT len);
static void idct_II(LC3_FLOAT* in, LC3_FLOAT* out, LC3_INT len);
static void pvq_subpyr_search(LC3_FLOAT* x_in, LC3_INT dim, LC3_INT pulses, LC3_INT* y, LC3_FLOAT* y_en1_norm);

void idct_II(LC3_FLOAT* in, LC3_FLOAT* out, LC3_INT len)
{
    LC3_INT   i = 0, j = 0;
    LC3_FLOAT norm1 = 0, norm2 = 0, sum = 0, sumInLoop = 0;

    norm1 = LC3_SQRT(2.0 / (LC3_FLOAT)len);
    norm2 = 1.0 / (LC3_SQRT(2.0));

    for (i = 0; i < len; i++) {
        sum = 0;
        for (j = 0; j < len; j++) {
            sumInLoop = in[j] * cos(M_PI / (2.0 * (LC3_FLOAT)len) * (2.0 * ((LC3_FLOAT)i + 1.0) - 1.0) * ((LC3_FLOAT)j));

            if (j == 0) {
                sumInLoop *= norm2;
            }

            sum += sumInLoop;
        }

        out[i] = norm1 * sum;
    }
}

void pvq_subpyr_search(LC3_FLOAT* x_in, LC3_INT dim, LC3_INT pulses, LC3_INT* y, LC3_FLOAT* y_en1_norm)
{
    LC3_FLOAT x[M] = {0}, xabs[M] = {0}, xsum = {0}, yy_tmp = 0, xy2_tmp = 0, gain_fac = 0;
    LC3_INT   i = 0, xsign[M] = {0}, imax = 0, pulse_tot = 0;
    LC3_FLOAT eps = LC3_POW(2, -24), proj_fac = 0, xy = 0, yy = 0, cmax_num = 0, cmax_den = 0;

    if (pulses == 0) {
        return;
    }

    move_float(x, x_in, dim);

    for (i = 0; i < dim; i++) {
        xabs[i] = fabs(x[i]);
    }

    for (i = 0; i < dim; i++) {
        if (x[i] >= 0) {
            xsign[i] = 1;
        } else {
            xsign[i] = -1;
        }
    }

    for (i = 0; i < dim; i++) {
        xsum += xabs[i];
    }

    if (xsum > eps) {
        pulse_tot = 0;
        yy        = 0;
        xy        = 0;

        /* Find a start position on a lower sub pyramid */
        proj_fac = (pulses - 1) / xsum;

        for (i = 0; i < dim; i++) {
            y[i] = floor(xabs[i] * proj_fac);
            pulse_tot += y[i];
            yy = yy + y[i] * y[i];
            xy = xy + xabs[i] * y[i];
        }

        /* Now run ACELP-like full corrsq/energy search */
        yy = yy * 0.5;

        while (pulse_tot < pulses) {
            imax     = 0;
            cmax_num = -LC3_POW(2, 15);
            cmax_den = 0;

            yy = yy + 0.5;

            for (i = 0; i < dim; i++) {
                xy2_tmp = xy + xabs[i];
                xy2_tmp = xy2_tmp * xy2_tmp;
                yy_tmp  = yy + y[i];

                if (xy2_tmp * cmax_den > yy_tmp * cmax_num) {
                    cmax_num = xy2_tmp;
                    cmax_den = yy_tmp;
                    imax     = i;
                }
            }

            xy        = xy + xabs[imax];
            yy        = yy + y[imax];
            y[imax]   = y[imax] + 1;
            pulse_tot = pulse_tot + 1;
        }

        yy = yy * 2.0;
    } else {
        pulse_tot = pulses;
        yy        = 0.0;

        if (dim > 1) {
            y[0]   = floor(pulses / 2);
            y[dim] = -(pulses - floor(pulses / 2));
            yy     = y[0] * y[0] + y[dim] * y[dim];
        } else {
            y[1] = pulses;
            yy   = pulses * pulses;
        }
    }

    /* Apply scaling to unit energy, always at least one pulse so no div-by-zero*/
    gain_fac = 1.0 * 1.0 / LC3_SQRT(yy);

    for (i = 0; i < dim; i++) {
        y[i]          = y[i] * xsign[i];
        y_en1_norm[i] = y[i] * gain_fac;
    }
}

void pvq_enc(LC3_INT* pulses, LC3_INT* LS_ind, LC3_INT* MPVQ_ind, LC3_INT len)
{
    LC3_INT k = 0, pos = 0;

    k         = 0;
    *LS_ind   = -1;
    *MPVQ_ind = 0;

    /* Encoding loop */
    for (pos = len - 1; pos >= 0; pos--) {
        if (*LS_ind >= 0 && pulses[pos] != 0) {
            *MPVQ_ind = 2 * (*MPVQ_ind) + *LS_ind;
        }

        if (pulses[pos] > 0) {
            *LS_ind = 0;
        }

        if (pulses[pos] < 0) {
            *LS_ind = 1;
        }

        *MPVQ_ind = *MPVQ_ind + pvq_enc_A[len - pos - 1][k];
        k         = k + abs(pulses[pos]);
    }
}

void process_snsQuantizesScf_Enc(LC3_FLOAT* env, LC3_INT* index, LC3_FLOAT* envq, Dct2 dct2structSNS)
{
    LC3_INT Ntot = 0, N1 = 0, sec = 0, st = 0, fin = 0, i = 0, j = 0, c = 0, idx = 0;

    LC3_FLOAT target[M] = {0}, sum = 0, bestdist = 0, st1_vector[M] = {0}, pvq_target_pre[M] = {0};
    LC3_FLOAT pvq_target[M]                = {0};
    LC3_FLOAT stage2_en1_norm_pre_subC1[M] = {0}, stage2_en1_norm_pre_sub[M] = {0}, tmp[M] = {0}, gain_fac = 0,
          yC_en1_norm[M] = {0}, normZero[M] = {0}, v[M][M] = {{0}}, min_err = 0, err_sig_split[M] = {0}, e1_split = 0;
    LC3_FLOAT pvq_glob_gain = 0, stage2_en1_norm_split[M] = {0}, st2_vector[M] = {0}, e1_sofar = 0,
          stage2_en1_norm_sub[M] = {0}, g = 0;

    LC3_INT pulses_submodeC1[M] = {0}, N2 = 0, Kh1 = 0, Kh2 = 0, pulses_submodeC2[M] = {0}, yC[M] = {0}, pulses[M] = {0},
        Ksub = 0, N_110 = 0, K_110 = 0, tmp1 = 0, tmp2 = 0;

    Ntot = 16;
    N1   = floor(Ntot / 2);

    /* Run Search */
    for (sec = 0; sec < 2; sec++) {
        st  = sec * N1 + 1;
        fin = st + N1 - 1;

        j = 0;
        for (i = st - 1; i < fin; i++) {
            target[j] = env[i];
            j++;
        }

        bestdist = LC3_POW(2, 100);

        for (c = 0; c < 32; c++) {
            if (sec == 0) {
                sum = 0;
                for (i = 0; i < (fin - st + 1); i++) {
                    sum += (target[i] - sns_C1[i][c]) * (target[i] - sns_C1[i][c]);
                }
            } else {
                sum = 0;
                for (i = 0; i < (fin - st + 1); i++) {
                    sum += (target[i] - sns_C2[i][c]) * (target[i] - sns_C2[i][c]);
                }
            }

            if (sum < bestdist) {
                bestdist = sum;
                idx      = c;
            }
        }

        index[sec] = idx;

        j = 0;
        for (i = st - 1; i < fin; i++) {
            if (sec == 0) {
                st1_vector[i] = sns_C1[j][idx];
            } else {
                st1_vector[i] = sns_C2[j][idx];
            }

            j++;
        }
    }

    /* STAGE 2 */

    for (i = 0; i < 16; i++) {
        pvq_target_pre[i] = env[i] - st1_vector[i];
    }

    dct2_apply(&dct2structSNS, pvq_target_pre, pvq_target);

    N1  = 10;
    N2  = Ntot - N1;
    Kh1 = 10;
    Kh2 = 1;

    pvq_subpyr_search(pvq_target, N1, Kh1, pulses_submodeC1, stage2_en1_norm_pre_subC1);

    pvq_subpyr_search(&pvq_target[10], N2, Kh2, pulses_submodeC2, tmp);

    move_int(yC, pulses_submodeC1, 10);

    move_int(&yC[10], pulses_submodeC2, 6);

    sum = 0;
    for (i = 0; i < M; i++) {
        sum += yC[i] * yC[i];
    }

    gain_fac = 1.0 / LC3_SQRT(sum);

    for (i = 0; i < M; i++) {
        yC_en1_norm[i] = yC[i] * gain_fac;
    }

    /* Gain */

    move_float(normZero, stage2_en1_norm_pre_subC1, 10);

    for (i = 0; i < M; i++) {
        v[0][i] = sns_vq_reg_adj_gains_fl[0] * yC_en1_norm[i];
        v[1][i] = sns_vq_reg_adj_gains_fl[1] * yC_en1_norm[i];
    }

    for (i = 0; i < M; i++) {
        v[2][i] = sns_vq_reg_lf_adj_gains_fl[0] * normZero[i];
        v[3][i] = sns_vq_reg_lf_adj_gains_fl[1] * normZero[i];
        v[4][i] = sns_vq_reg_lf_adj_gains_fl[2] * normZero[i];
        v[5][i] = sns_vq_reg_lf_adj_gains_fl[3] * normZero[i];
    }

    min_err = LC3_POW(2, 15);

    for (i = 0; i < 6; i++) {
        sum = 0;
        for (j = 0; j < M; j++) {
            sum += (pvq_target[j] - v[i][j]) * (pvq_target[j] - v[i][j]);
        }

        if (sum < min_err) {
            min_err = sum;
            idx     = i;
        }
    }

    for (i = 0; i < M; i++) {
        yC_en1_norm[i] = v[idx][i] / q_g_sns[idx];
    }

    pvq_glob_gain = q_g_sns[idx];

    /* */
    idct_II(yC_en1_norm, stage2_en1_norm_split, M);

    /* Distortion */
    sum = 0;
    for (i = 0; i < M; i++) {
        err_sig_split[i] = pvq_target_pre[i] - pvq_glob_gain * stage2_en1_norm_split[i];
        sum += err_sig_split[i] * err_sig_split[i];
    }

    e1_split = sum;

    /* Mode selection */
    e1_sofar = LC3_POW(2, 15);

    if (e1_split < e1_sofar) {
        if (idx <= 1) {
            index[2] = 0;
            index[3] = idx;

            move_int(pulses, yC, M);

        } else {
            index[2] = 1;
            index[3] = idx - 2;

            move_int(pulses, pulses_submodeC1, 10);
        }

        for (i = 0; i < M; i++) {
            st2_vector[i] = pvq_glob_gain * stage2_en1_norm_split[i];
        }

        e1_sofar = e1_split;
    }

    /* Outlier near mode */
    Ksub = 8;
    pvq_subpyr_search(pvq_target, Ntot, Ksub, yC, stage2_en1_norm_pre_sub);

    idct_II(stage2_en1_norm_pre_sub, stage2_en1_norm_sub, M);

    /* Gain */
    min_err = LC3_POW(2, 15);

    for (i = 0; i < 4; i++) {
        g = sns_vq_near_adj_gains_fl[i];

        sum = 0;
        for (j = 0; j < M; j++) {
            sum += (pvq_target_pre[j] - g * stage2_en1_norm_sub[j]) * (pvq_target_pre[j] - g * stage2_en1_norm_sub[j]);
        }

        if (sum < min_err) {
            idx           = i;
            min_err       = sum;
            pvq_glob_gain = g;
        }
    }

    /* Mode selection */
    if (min_err < e1_sofar) {
        index[2] = 2;
        index[3] = idx;

        for (i = 0; i < M; i++) {
            st2_vector[i] = pvq_glob_gain * stage2_en1_norm_sub[i];
        }

        move_int(pulses, yC, M);

        e1_sofar = min_err;
    }

    /* Outlier far mode */
    N_110 = Ntot;
    K_110 = 6;

    /* Pulse search */
    pvq_subpyr_search(pvq_target, N_110, K_110, yC, stage2_en1_norm_pre_sub);

    /* Inverse transform */
    idct_II(stage2_en1_norm_pre_sub, stage2_en1_norm_sub, M);

    /* Gain */
    min_err = LC3_POW(2, 15);

    for (i = 0; i < 8; i++) {
        g = sns_vq_far_adj_gains_fl[i];

        sum = 0;
        for (j = 0; j < M; j++) {
            sum += (pvq_target_pre[j] - g * stage2_en1_norm_sub[j]) * (pvq_target_pre[j] - g * stage2_en1_norm_sub[j]);
        }

        if (sum < min_err) {
            idx           = i;
            min_err       = sum;
            pvq_glob_gain = g;
        }
    }

    /* Mode selection */

    if (min_err < e1_sofar) {
        index[2] = 3;
        index[3] = idx;

        for (i = 0; i < M; i++) {
            pulses[i]     = yC[i];
            st2_vector[i] = pvq_glob_gain * stage2_en1_norm_sub[i];
        }
    }

    /* Indexing */

    if (index[2] < 2) {
        pvq_enc(pulses, &index[4], &index[5], 10);
    } else {
        pvq_enc(pulses, &index[4], &index[5], M);
    }

    if (index[2] == 0) {
        pvq_enc(&pulses[10], &tmp1, &tmp2, 6);
        index[6] = tmp2 * 2 + tmp1;
    } else if (index[2] == 2) {
        index[6] = -1;
    } else {
        index[6] = -2;
    }

    for (i = 0; i < M; i++) {
        envq[i] = st1_vector[i] + st2_vector[i];
    }
}

LC3_INT find_last_indice_le(LC3_INT compare, const LC3_INT* array, LC3_INT len)
{
    LC3_INT idx = 0, i = 0;

    for (i = 0; i < len; i++) {
        if (compare >= array[i]) {
            idx++;
        }
    }

    if (idx > 0) {
        idx--;
    }

    return idx;
}

void pvq_dec(LC3_INT k, LC3_INT m, LC3_INT LS_ind, LC3_INT MPVQ_ind, LC3_INT* pulses)
{
    LC3_INT leading_sign = 0, idx = 0, k_delta = 0, pos = 0;

    leading_sign = 1 - 2 * LS_ind;

    /* Decoding loop */

    for (pos = 0; pos < m; pos++) {
        if (MPVQ_ind != 0) {
            /* Find last indice */
            idx      = find_last_indice_le(MPVQ_ind, &pvq_enc_A[m - pos - 1][0], k + 1);
            MPVQ_ind = MPVQ_ind - pvq_enc_A[m - pos - 1][idx];
            k_delta  = k - idx;
        } else {
            pulses[pos] = leading_sign * k;
            break;
        }

        if (k_delta != 0) {
            pulses[pos] = leading_sign * k_delta;
            if ((MPVQ_ind % 2) != 0) {
                leading_sign = -1;
            } else {
                leading_sign = 1;
            }

            MPVQ_ind = floor(MPVQ_ind / 2);
            k        = k - k_delta;
        }
    }
}

void process_snsQuantizesScf_Dec(LC3_INT* scf_idx, LC3_FLOAT* scf_q)
{
    LC3_INT   i = 0, submode = 0;
    LC3_INT   pulses2[6] = {0}, pulses[M] = {0};
    LC3_FLOAT st2_vector[M] = {0}, st2_vector_idct[M] = {0}, sum = 0;

    /* Decode first stage */

    for (i = 0; i < 8; i++) {
        scf_q[i]     = sns_C1[i][scf_idx[0]];
        scf_q[i + 8] = sns_C2[i][scf_idx[1]];
    }

    /* STAGE 2 */
    /* Decode submode */

    submode = scf_idx[2];

    /* Decode pulses */

    if (submode < 2) {
        pvq_dec(10, 10, scf_idx[4], scf_idx[5], pulses);

        if (submode == 0) {
            pvq_dec(1, 6, (scf_idx[6] % 2), floor(scf_idx[6] / 2), pulses2);

            move_int(&pulses[10], pulses2, 6);

        } else {
            pulses[15] = 0;
        }
    } else if (submode == 2) {
        pvq_dec(8, 16, scf_idx[4], scf_idx[5], pulses);
    } else {
        pvq_dec(6, 16, scf_idx[4], scf_idx[5], pulses);
    }

    /* Normalization */

    for (i = 0; i < M; i++) {
        sum += pulses[i] * pulses[i];
    }

    sum = LC3_SQRT(sum);

    for (i = 0; i < M; i++) {
        st2_vector[i] = pulses[i] / sum;
    }

    /* Inverse transform */
    idct_II(st2_vector, st2_vector_idct, M);

    /* Gain */
    for (i = 0; i < M; i++) {
        st2_vector_idct[i] = st2_vector_idct[i] * sns_dec_gains[submode][scf_idx[3]];
    }

    /* Add stage 1 and stage 2 */

    for (i = 0; i < M; i++) {
        scf_q[i] = scf_q[i] + st2_vector_idct[i];
    }
}
