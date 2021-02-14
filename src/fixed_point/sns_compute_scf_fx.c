/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"



void processSnsComputeScf_fx(Word32 *d2_fx, Word16 d2_fx_exp, Word16 fs_idx, Word16 n_bands, Word16 *scf,
                             Word16 scf_smoothing_enabled, Word8 *scratchBuffer)
{
    Dyn_Mem_Deluxe_In(
        Word16  i, s, s2, nf;
        Word32  L_mean, L_tmp;
        Word32 *d3_fx;
        Word16 *d3_fx_exp;
        Word16 *d4_fx;
        Word16 *scf_smooth;
    );


    d3_fx     = scratchAlign(scratchBuffer, 0);                         /* Size = 4 * MAX_BANDS_NUMBER = 256 bytes */
    d3_fx_exp = scratchAlign(d3_fx, sizeof(*d3_fx) * MAX_BANDS_NUMBER); /* Size = 2 * MAX_BANDS_NUMBER = 128 bytes */
    d4_fx = scratchAlign(d3_fx_exp, sizeof(*d3_fx_exp) * MAX_BANDS_NUMBER); /* Size = 2 * MAX_BANDS_NUMBER = 128bytes */
    scf_smooth = scratchAlign(d4_fx, sizeof(*d4_fx) * MAX_BANDS_NUMBER);    /* Size = 2 * 16 */

/* Smoothing and Pre-emphasis */
    IF (sub(n_bands, 32) < 0)
    {
        L_tmp = sub(32, n_bands);
        FOR (i = sub(n_bands, 1); i >= L_tmp; i--)
        {
            d2_fx[(i + L_tmp) * 2 + 1] = d2_fx[i]; move32();
            d2_fx[(i + L_tmp) * 2 + 0] = d2_fx[i]; move32();
        }
        FOR (i = sub(L_tmp, 1); i >= 0; i--)
        {
            d2_fx[i * 4 + 3] = d2_fx[i]; move32();
            d2_fx[i * 4 + 2] = d2_fx[i]; move32();
            d2_fx[i * 4 + 1] = d2_fx[i]; move32();
            d2_fx[i * 4 + 0] = d2_fx[i]; move32();
        }
        n_bands = 64; move16();
    }
    ELSE
    IF (sub(n_bands, 64) < 0)
    {
        L_tmp = sub(64, n_bands);
        FOR (i = sub(n_bands, 1); i >= L_tmp; i--)
        {
            d2_fx[i + L_tmp] = d2_fx[i]; move32();
        }
        FOR (i = sub(L_tmp, 1); i >= 0; i--)
        {
            d2_fx[i * 2 + 1] = d2_fx[i]; move32();
            d2_fx[i * 2 + 0] = d2_fx[i]; move32();
        }
        n_bands = 64; move16();
    }
    
    L_tmp        = L_add(Mpy_32_16(d2_fx[0], 24576), L_shr_pos(d2_fx[1], 2));
    d3_fx[0]     = Mpy_32_16(L_tmp, lpc_pre_emphasis[fs_idx][0]); move32();
    d3_fx_exp[0] = add(d2_fx_exp, lpc_pre_emphasis_e[fs_idx][0]); move16();
    FOR (i = 1; i < n_bands - 1; i++)
    {
        L_tmp        = L_add(L_shr_pos(d2_fx[i], 1), L_add(L_shr_pos(d2_fx[i - 1], 2), L_shr_pos(d2_fx[i + 1], 2)));
        d3_fx[i]     = Mpy_32_16(L_tmp, lpc_pre_emphasis[fs_idx][i]); move32();
        d3_fx_exp[i] = add(d2_fx_exp, lpc_pre_emphasis_e[fs_idx][i]); move16();
    }
    L_tmp                  = L_add(Mpy_32_16(d2_fx[n_bands - 1], 24576), L_shr_pos(d2_fx[n_bands - 2], 2));
    d3_fx[n_bands - 1]     = Mpy_32_16(L_tmp, lpc_pre_emphasis[fs_idx][n_bands - 1]); move32();
    d3_fx_exp[n_bands - 1] = add(d2_fx_exp, lpc_pre_emphasis_e[fs_idx][n_bands - 1]); move16();

    /* Mean */
    s  = d3_fx_exp[MAX_BANDS_NUMBER - 1];
    s2 = add(s, 6);

    L_mean = L_shr(d3_fx[0], sub(s2, d3_fx_exp[0]));
    FOR (i = 1; i < MAX_BANDS_NUMBER; i++)
    {
        L_mean = L_add(L_mean, L_shr(d3_fx[i], sub(s2, d3_fx_exp[i])));
    }

    /* Noise floor at -40dB */
    nf = BASOP_Util_Log2_16(L_mean, s);
    nf = sub(s_max(nf, -25965), 6803);

/* Log-domain */
    FOR (i = 0; i < MAX_BANDS_NUMBER; i++)
    {
        d4_fx[i] = s_max(nf, BASOP_Util_Log2_16(d3_fx[i], d3_fx_exp[i])); move16();
    }

    /* Downsampling */
    L_tmp    = L_mult(d4_fx[0], 8192);
    L_tmp    = L_mac(L_tmp, d4_fx[1], 8192);
    L_tmp    = L_mac(L_tmp, d4_fx[2], 8192);
    L_tmp    = L_mac(L_tmp, d4_fx[3], 5461);
    d3_fx[0] = L_mac(L_tmp, d4_fx[4], 2731); move32();
    FOR (i = 1; i < M - 1; i++)
    {
        L_tmp    = L_mult(d4_fx[i * 4 - 1], 2731);
        L_tmp    = L_mac(L_tmp, d4_fx[i * 4 + 0], 5461);
        L_tmp    = L_mac(L_tmp, d4_fx[i * 4 + 1], 8192);
        L_tmp    = L_mac(L_tmp, d4_fx[i * 4 + 2], 8192);
        L_tmp    = L_mac(L_tmp, d4_fx[i * 4 + 3], 5461);
        d3_fx[i] = L_mac(L_tmp, d4_fx[i * 4 + 4], 2731); move32();
    }
    L_tmp        = L_mult(d4_fx[59], 2731);
    L_tmp        = L_mac(L_tmp, d4_fx[60], 5461);
    L_tmp        = L_mac(L_tmp, d4_fx[61], 8192);
    L_tmp        = L_mac(L_tmp, d4_fx[62], 8192);
    d3_fx[M - 1] = L_mac(L_tmp, d4_fx[63], 8192); move32();

/* Remove mean and scaling */
    L_mean = L_shr_pos(d3_fx[0], 4);
    FOR (i = 1; i < M; i++)
    {
        L_mean = L_add(L_mean, L_shr_pos(d3_fx[i], 4));
    }

    FOR (i = 0; i < M; i++)
    {
        scf[i] = mult_r(27853, round_fx(L_shl_pos(L_sub(d3_fx[i], L_mean), 1))); move16();
    }

    /* scale factor smoothing */
    IF (scf_smoothing_enabled)
    {
        scf_smooth[0] = L_shr(L_mult0(L_add(L_add(scf[0], scf[1]), scf[2]), 10923), 15);
        L_mean        = scf_smooth[0]; move16();
        scf_smooth[1] = L_shr(L_add(L_add(L_add(scf[0], scf[1]), scf[2]), scf[3]), 2);
        L_mean        = L_add(L_mean, scf_smooth[1]);
        FOR (i = 2; i < M - 2; i++)
        {
            L_tmp         = L_add(L_add(L_add(L_add(scf[i - 2], scf[i - 1]), scf[i]), scf[i + 1]), scf[i + 2]);
            scf_smooth[i] = L_shr(L_mult0(L_tmp, 13107), 16);
            L_mean        = L_add(L_mean, scf_smooth[i]);
        }
        scf_smooth[M - 2] = L_shr(L_add(L_add(L_add(scf[M - 4], scf[M - 3]), scf[M - 2]), scf[M - 1]), 2);
        L_mean            = L_add(L_mean, scf_smooth[M - 2]);
        scf_smooth[M - 1] = L_shr(L_mult0(L_add(L_add(scf[M - 3], scf[M - 2]), scf[M - 1]), 10923), 15);
        L_mean            = L_add(L_mean, scf_smooth[M - 1]);

        L_mean = L_shr(L_mean, 4); // assumes M = 16
        FOR (i = 0; i < M; i++)
        {
            scf[i] = L_shr(L_sub(scf_smooth[i], L_mean), 1);
        }
    }

    Dyn_Mem_Deluxe_Out();
}

