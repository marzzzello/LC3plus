/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"



void processSnsInterpolateScf_fx(Word16 *scf_q, Word16 mdct_scf[], Word16 mdct_scf_exp[], Word16 inv_scf,
                                 Word16 n_bands, Word8 *scratchBuffer)
{
    Dyn_Mem_Deluxe_In(
        Word16  i, tmp2;
        Word16 *scf_int;
        Word16  tmp;
        Word16 *scf_tmp;
    );

    scf_int = (Word16 *)scratchAlign(scratchBuffer, 0);              /* Size = 2 * MAX_BANDS_NUMBER = 128 bytes */
    scf_tmp = (Word16 *)scratchAlign(scf_int, 2 * MAX_BANDS_NUMBER); /* 2 * MAX_BANDS_NUMBER = 128 bytes */

    /* Interpolation */
    scf_int[0] = scf_q[0];
    scf_int[1] = scf_q[0];
    FOR (i = 1; i < M; i++)
    {
        tmp                = sub(scf_q[i], scf_q[i - 1]);
        tmp2               = mult_r(tmp, 8192);
        tmp                = mult_r(tmp, 4096);
        scf_int[i * 4 - 2] = add(scf_q[i - 1], tmp);
        scf_int[i * 4 - 1] = add(scf_int[i * 4 - 2], tmp2);
        scf_int[i * 4]     = add(scf_int[i * 4 - 1], tmp2);
        scf_int[i * 4 + 1] = add(scf_int[i * 4], tmp2);
    }
    scf_int[62] = add(scf_int[61], tmp2);
    scf_int[63] = add(scf_int[62], tmp2);

    /* 8 kHz mode for 2.5 ms */
    IF (sub(n_bands, 32) < 0)
    {
        basop_memmove(scf_tmp, scf_int, 64 * sizeof(Word16));
        tmp = sub(32, n_bands);
        FOR (i = 0; i < tmp; i++)
        {
            /* 8192 = 0.25 * 2^15 */
            scf_int[i] = add(mac_r(L_mult(scf_tmp[4 * i], 8192), scf_tmp[4 * i + 1], 8192),
                             mac_r(L_mult(scf_tmp[4 * i + 2], 8192), scf_tmp[4 * i + 3], 8192));
        }

        FOR (i = 0; i < n_bands - tmp; i++)
        {
            scf_int[tmp + i] = mac_r(L_mult(scf_tmp[4 * tmp + 2 * i], 16384), scf_tmp[4 * tmp + 2 * i + 1], 16384);
        }
    }
    ELSE
        /* For 5ms */
        IF (sub(n_bands, 64) < 0)
    {
        tmp = sub(64, n_bands);
        FOR (i = 0; i < tmp; i++)
        {
            scf_int[i] = mac_r(L_mult(scf_int[2 * i], 16384), scf_int[2 * i + 1], 16384);
        }
        FOR (; i < n_bands; i++)
        {
            scf_int[i] = scf_int[tmp + i];
        }
    }

    /* Inversion at encoder-side*/
    IF (inv_scf == 1)
    {
        FOR (i = 0; i < n_bands; i++)
        {
            scf_int[i] = negate(scf_int[i]);
        }
    }

/* Linear Domain */
    FOR (i = 0; i < n_bands; i++)
    {
        mdct_scf[i] = BASOP_Util_InvLog2_16(scf_int[i], &mdct_scf_exp[i]);
    }

    Dyn_Mem_Deluxe_Out();
}

