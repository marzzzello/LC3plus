/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void process_ltpf_decoder_fl(LC3_FLOAT* x, LC3_INT xLen, LC3_FLOAT* y, LC3_INT fs, LC3_FLOAT* mem_old_x, LC3_FLOAT* mem_old_y,
                             LC3_INT* mem_pitch_int, LC3_INT* mem_pitch_fr, LC3_FLOAT* mem_gain, LC3_INT* mem_beta_idx, LC3_INT bfi,
                             LC3_INT* param, LC3_INT* mem_param, LC3_INT conf_beta_idx, LC3_FLOAT conf_beta, 
                             LC3_INT concealMethod, 
                             LC3_FLOAT damping)
{
    LC3_INT i = 0, j = 0, n = 0, conf_pitmin = 0, conf_pitfr2 = 0, conf_pitfr1 = 0, N = 0, L_past_x = 0, N4 = 0, N34 = 0,
        pitch_int = 0, pitch_fr = 0, prev_param[3] = {0}, p1 = 0, p2 = 0, L_past_y = 0, inter_len = 0, tilt_len = 0,
        tilt_len_r = 0, inter_len_r = 0, old_x_len = 0, old_y_len = 0;

    LC3_FLOAT conf_alpha = 0, gain = 0, a1[MAX_LEN] = {0}, a2[MAX_LEN] = {0}, b1[MAX_LEN] = {0}, b2[MAX_LEN] = {0},
          buf_x[4 * MAX_LEN] = {0}, buf_y[4 * MAX_LEN] = {0}, buf_z[4 * MAX_LEN] = {0}, pitch = 0, sum1 = 0, sum2 = 0;

    const LC3_FLOAT *inter_filter[4], *tilt_filter[4];

    conf_pitmin = 32;
    conf_pitfr2 = 127;
    conf_pitfr1 = 157;
    conf_alpha  = 0.85;

    if (fs == 8000 || fs == 16000) {
        inter_filter[0] = conf_inter_filter_16[0];
        inter_filter[1] = conf_inter_filter_16[1];
        inter_filter[2] = conf_inter_filter_16[2];
        inter_filter[3] = conf_inter_filter_16[3];
        inter_len_r     = 4;

        tilt_filter[0] = conf_tilt_filter_16[0];
        tilt_filter[1] = conf_tilt_filter_16[1];
        tilt_filter[2] = conf_tilt_filter_16[2];
        tilt_filter[3] = conf_tilt_filter_16[3];
        tilt_len       = 4 - 2;
        tilt_len_r     = 3;
    } else if (fs == 24000) {
        inter_filter[0] = conf_inter_filter_24[0];
        inter_filter[1] = conf_inter_filter_24[1];
        inter_filter[2] = conf_inter_filter_24[2];
        inter_filter[3] = conf_inter_filter_24[3];
        inter_len_r     = 6;

        tilt_filter[0] = conf_tilt_filter_24[0];
        tilt_filter[1] = conf_tilt_filter_24[1];
        tilt_filter[2] = conf_tilt_filter_24[2];
        tilt_filter[3] = conf_tilt_filter_24[3];
        tilt_len       = 6 - 2;
        tilt_len_r     = 5;
    } else if (fs == 32000) {
        inter_filter[0] = conf_inter_filter_32[0];
        inter_filter[1] = conf_inter_filter_32[1];
        inter_filter[2] = conf_inter_filter_32[2];
        inter_filter[3] = conf_inter_filter_32[3];
        inter_len_r     = 8;

        tilt_filter[0] = conf_tilt_filter_32[0];
        tilt_filter[1] = conf_tilt_filter_32[1];
        tilt_filter[2] = conf_tilt_filter_32[2];
        tilt_filter[3] = conf_tilt_filter_32[3];
        tilt_len       = 8 - 2;
        tilt_len_r     = 7;
    } else if (fs == 44100 || fs == 48000) {
        inter_filter[0] = conf_inter_filter_48[0];
        inter_filter[1] = conf_inter_filter_48[1];
        inter_filter[2] = conf_inter_filter_48[2];
        inter_filter[3] = conf_inter_filter_48[3];
        inter_len_r     = 12;

        tilt_filter[0] = conf_tilt_filter_48[0];
        tilt_filter[1] = conf_tilt_filter_48[1];
        tilt_filter[2] = conf_tilt_filter_48[2];
        tilt_filter[3] = conf_tilt_filter_48[3];
        tilt_len       = 12 - 2;
        tilt_len_r     = 11;
    }

    inter_len = MAX(fs, 16000) / 8000;

    /* Init buffers */
    N         = xLen;
    old_x_len = tilt_len;
    old_y_len = ceil(228 * fs / 12800) + inter_len;

    L_past_x = old_x_len;

    move_float(buf_x, mem_old_x, old_x_len);

    move_float(&buf_x[old_x_len], x, xLen);

    L_past_y = old_y_len;

    move_float(buf_y, mem_old_y, old_y_len);

    zero_float(&buf_y[old_y_len], xLen);

    N4  = N / 4;
    N34 = 3 * N / 4;

    if (bfi && concealMethod == 0) {
        param[0] = 0;
        param[1] = 0;
    }

    if (!bfi) {
        /* Decode pitch */
        if (param[0] == 1) {
            if (param[2] < (conf_pitfr2 - conf_pitmin) * 4) {
                pitch_int = conf_pitmin + floor(param[2] / 4);
                pitch_fr  = param[2] - ((pitch_int - conf_pitmin) * 4);
            } else if (param[2] < ((conf_pitfr2 - conf_pitmin) * 4) + ((conf_pitfr1 - conf_pitfr2) * 2)) {
                param[2]  = param[2] - ((conf_pitfr2 - conf_pitmin) * 4);
                pitch_int = conf_pitfr2 + floor(param[2] / 2);
                pitch_fr  = param[2] - ((pitch_int - conf_pitfr2) * 2);
                pitch_fr  = pitch_fr * 2;
            } else {
                pitch_int =
                    param[2] + (conf_pitfr1 - ((conf_pitfr2 - conf_pitmin) * 4) - ((conf_pitfr1 - conf_pitfr2) * 2));
                pitch_fr = 0;
            }

            pitch     = ((LC3_FLOAT)pitch_int + (LC3_FLOAT)pitch_fr / 4.0) * (LC3_FLOAT)fs / 12800.0;
            pitch     = round(pitch * 4.0) / 4.0;
            pitch_int = floor(pitch);
            pitch_fr  = (LC3_INT)((pitch - (LC3_FLOAT)pitch_int) * 4.0);
        } else {
            pitch_int = 0;
            pitch_fr  = 0;
        }

        /* Decode gain */
        if (conf_beta_idx < 0) {
            param[1] = 0;
        }

        if (param[1] == 1) {
            gain = conf_beta;
        } else {
            gain = 0;
        }
    } else {
        if (conf_beta_idx < 0) {
            if (mem_param[1] && *mem_beta_idx >= 0)
            {
                conf_beta_idx = *mem_beta_idx;
            }
        }
        
        memmove(param, mem_param, sizeof(LC3_INT) * 3);
        
        if (concealMethod == 2)
        {
            /* cause the ltpf to "fade_out" and only filter during initial 2.5 ms and then its buffer during 7.5 ms */
            assert(bfi == 1);
            param[1] = 0; /* ltpf_active = 0 */
        }
        
        pitch_int = *mem_pitch_int;
        pitch_fr  = *mem_pitch_fr;
        gain      = (LC3_FLOAT) *mem_gain * damping;
    }

    /* Init filter parameters */

    memmove(prev_param, mem_param, sizeof(LC3_INT) * 3);

    if (prev_param[1] == 1) {
        for (i = 0; i < inter_len_r; i++) {
            a1[i] = *mem_gain * inter_filter[*mem_pitch_fr][i];
        }

        for (i = 0; i < tilt_len_r; i++) {
            b1[i] = conf_alpha * (*mem_gain) * tilt_filter[*mem_beta_idx][i];
        }

        p1 = *mem_pitch_int;
    }

    if (param[1] == 1) {
        for (i = 0; i < tilt_len_r; i++) {
            b2[i] = conf_alpha * gain * tilt_filter[conf_beta_idx][i];
        }

        for (i = 0; i < inter_len_r; i++) {
            a2[i] = gain * inter_filter[pitch_fr][i];
        }

        p2 = pitch_int;
    }

    /* First quarter of the current frame: cross-fading */
    if (prev_param[1] == 0 && param[1] == 0) {
        memmove(&buf_y[L_past_y], &buf_x[L_past_x], sizeof(LC3_FLOAT) * N4);

    } else if (prev_param[1] == 1 && param[1] == 0) {
        for (n = 0; n < N4; n++) {
            sum1 = 0;
            sum2 = 0;
            j    = 0;
            for (i = L_past_x + n; i >= L_past_x + n - tilt_len; i--) {
                sum1 += b1[j] * buf_x[i];
                j++;
            }

            j = 0;
            for (i = L_past_y + n - p1 + inter_len - 1; i >= L_past_y + n - p1 - inter_len; i--) {
                sum2 += a1[j] * buf_y[i];
                j++;
            }

            buf_y[L_past_y + n] = buf_x[L_past_x + n] - (((LC3_FLOAT)N4 - (LC3_FLOAT)n) / (LC3_FLOAT)N4) * sum1 +
                                  (((LC3_FLOAT)N4 - (LC3_FLOAT)n) / (LC3_FLOAT)N4) * sum2;
        }

    } else if (prev_param[1] == 0 && param[1] == 1) {
        for (n = 0; n < N4; n++) {
            sum1 = 0;
            sum2 = 0;
            j    = 0;
            for (i = L_past_x + n; i >= L_past_x + n - tilt_len; i--) {
                sum1 += b2[j] * buf_x[i];
                j++;
            }

            j = 0;
            for (i = L_past_y + n - p2 + inter_len - 1; i >= L_past_y + n - p2 - inter_len; i--) {
                sum2 += a2[j] * buf_y[i];
                j++;
            }

            buf_y[L_past_y + n] = buf_x[L_past_x + n] - ((LC3_FLOAT)n / (LC3_FLOAT)N4) * sum1 + ((LC3_FLOAT)n / (LC3_FLOAT)N4) * sum2;
        }
    } else if (*mem_pitch_int == pitch_int && *mem_pitch_fr == pitch_fr) {
        for (n = 0; n < N4; n++) {
            sum1 = 0;
            sum2 = 0;
            j    = 0;
            for (i = L_past_x + n; i >= L_past_x + n - tilt_len; i--) {
                sum1 += b2[j] * buf_x[i];
                j++;
            }

            j = 0;
            for (i = L_past_y + n - p2 + inter_len - 1; i >= L_past_y + n - p2 - inter_len; i--) {
                sum2 += a2[j] * buf_y[i];
                j++;
            }

            buf_y[L_past_y + n] = buf_x[L_past_x + n] - sum1 + sum2;
        }
    } else {
        for (n = 0; n < N4; n++) {
            sum1 = 0;
            sum2 = 0;
            j    = 0;
            for (i = L_past_x + n; i >= L_past_x + n - tilt_len; i--) {
                sum1 += b1[j] * buf_x[i];
                j++;
            }

            j = 0;
            for (i = L_past_y + n - p1 + inter_len - 1; i >= L_past_y + n - p1 - inter_len; i--) {
                sum2 += a1[j] * buf_y[i];
                j++;
            }

            buf_y[L_past_y + n] = buf_x[L_past_x + n] - (((LC3_FLOAT)N4 - (LC3_FLOAT)n) / (LC3_FLOAT)N4) * sum1 +
                                  (((LC3_FLOAT)N4 - (LC3_FLOAT)n) / (LC3_FLOAT)N4) * sum2;
        }

        memmove(buf_z, buf_y, sizeof(LC3_FLOAT) * 4 * MAX_LEN);

        for (n = 0; n < N4; n++) {
            sum1 = 0;
            sum2 = 0;
            j    = 0;
            for (i = L_past_y + n; i >= L_past_y + n - tilt_len; i--) {
                sum1 += b2[j] * buf_z[i];
                j++;
            }

            j = 0;
            for (i = L_past_y + n - p2 + inter_len - 1; i >= L_past_y + n - p2 - inter_len; i--) {
                sum2 += a2[j] * buf_y[i];
                j++;
            }

            buf_y[L_past_y + n] = buf_z[L_past_y + n] - ((LC3_FLOAT)n / (LC3_FLOAT)N4) * sum1 + ((LC3_FLOAT)n / (LC3_FLOAT)N4) * sum2;
        }
    }

    /* Second quarter of the current frame */
    if (param[1] == 0) {
        memmove(&buf_y[L_past_y + N4], &buf_x[L_past_x + N4],
                sizeof(LC3_FLOAT) * ((L_past_x + N4 + N34) - (L_past_x + N4)));
    } else {
        for (n = 0; n < N34; n++) {
            sum1 = 0;
            sum2 = 0;
            j    = 0;
            for (i = L_past_x + N4 + n; i >= L_past_x + n + N4 - tilt_len; i--) {
                sum1 += b2[j] * buf_x[i];
                j++;
            }

            j = 0;
            for (i = L_past_y + N4 + n - p2 + inter_len - 1; i >= L_past_y + N4 + n - p2 - inter_len; i--) {
                sum2 += a2[j] * buf_y[i];
                j++;
            }

            buf_y[L_past_y + N4 + n] = buf_x[L_past_x + N4 + n] - sum1 + sum2;
        }
    }

    /* Update memory */

    move_float(mem_old_x, &buf_x[N], old_x_len);

    move_float(mem_old_y, &buf_y[N], old_y_len);

    move_int(mem_param, param, 3);

    *mem_pitch_int = pitch_int;
    *mem_pitch_fr  = pitch_fr;
    *mem_gain      = gain;
    *mem_beta_idx  = conf_beta_idx;

    /* Output */

    move_float(y, &buf_y[L_past_y], N);
}
