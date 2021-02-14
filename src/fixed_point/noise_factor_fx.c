/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"


#ifdef NONBE_LOW_BR_NF_TUNING
void processNoiseFactor_fx(Word16 *fac_ns_idx, Word16 x_e, Word32 x[], Word16 xq[], Word16 gg, Word16 gg_e,
                           Word16 BW_cutoff_idx, Word16 frame_dms, Word16 target_bytes, Word8 *scratchBuffer)
#else
void processNoiseFactor_fx(Word16 *fac_ns_idx, Word16 x_e, Word32 x[], Word16 xq[], Word16 gg, Word16 gg_e,
                           Word16 BW_cutoff_idx, Word16 frame_dms, Word8 *scratchBuffer)
#endif
{
    Dyn_Mem_Deluxe_In(
        Counter k;
        Word16  nzeros, s1, s2, s3, c, idx, fac_unq, *ind;
        Word16  noisefillwidth, noisefillstart, N;
        Word32  Lsum;
    );

    ind = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * MAX_LEN bytes */

    noisefillwidth = 0;
    noisefillstart = 0;
    c              = 0;                                move16();
    N              = BW_cutoff_bin_all[BW_cutoff_idx]; move16();

    SWITCH (frame_dms)
    {
    case 25:
        N              = shr_pos(N, 2);
        noisefillwidth = NOISEFILLWIDTH_2_5MS;
        noisefillstart = NOISEFILLSTART_2_5MS;
        BREAK;
    case 50:
        N              = shr_pos(N, 1);
        noisefillwidth = NOISEFILLWIDTH_5MS;
        noisefillstart = NOISEFILLSTART_5MS;
        BREAK;
    case 100:
        noisefillwidth = NOISEFILLWIDTH;
        noisefillstart = NOISEFILLSTART;
        BREAK;
    }

    nzeros = -2 * noisefillwidth - 1; move16();

    FOR (k = noisefillstart - noisefillwidth; k < noisefillstart + noisefillwidth; k++)
    {
        if (xq[k] != 0)
        {
            nzeros = -2 * noisefillwidth - 1; move16();
        }
        if (xq[k] == 0)
        {
            nzeros = add(nzeros, 1);
        }
    }

    FOR (k = noisefillstart; k < N - noisefillwidth; k++)
    {
        if (xq[k + noisefillwidth] != 0)
        {
            nzeros = -2 * noisefillwidth - 1; move16();
        }
        if (xq[k + noisefillwidth] == 0)
        {
            nzeros = add(nzeros, 1);
        }
        if (nzeros >= 0)
        {
            ind[c++] = k; move16();
        }
    }

    FOR (k = N - noisefillwidth; k < N; k++)
    {
        nzeros = add(nzeros, 1);
        if (nzeros >= 0)
        {
            ind[c++] = k; move16();
        }
    }

    IF (c == 0)
    {
        fac_unq = 0; move16();
    }
    ELSE
    {

#ifdef NONBE_LOW_BR_NF_TUNING
        IF (target_bytes <= 20 && frame_dms == 100)
        {
            Word32 ind_sum;
            Word16 mean_ind;

            Word16 fac_unq1, fac_unq2;

            /* calculate mean index */
            ind_sum = ind[0]; move32();
            FOR (k = 1; k < c; k++)
            {
                ind_sum = L_add(ind_sum, ind[k]);
            }

            mean_ind = BASOP_Util_Divide3216_Scale(ind_sum, c, &s2);
            mean_ind = shl(mean_ind, s2 + 1);

            assert(0 <= mean_ind && mean_ind <= ind[c - 1]);

            /* calculate noise filling gain for low frequencies */
            Lsum = 0; move32();
            FOR (k = 0; ind[k] <= mean_ind; k++)
            {
                Lsum = L_add(Lsum, L_abs(x[ind[k]]));
            }
            fac_unq1 = BASOP_Util_Divide3216_Scale(Lsum, k, &s1);
            fac_unq1 = BASOP_Util_Divide1616_Scale(fac_unq1, gg, &s2);
            s3       = sub(15, add(x_e, add(s1, sub(s2, gg_e))));
            s2       = norm_s(fac_unq1);
            test();
            IF (fac_unq1 != 0 && add(s3, s2) < 0)
            {
                fac_unq1 = MAX_16; move16();
            }
            ELSE
            {
                fac_unq1 = shr_r(fac_unq1, s_min(s3, 15));
            }

            /* calculate noise filling gain for high frequencies */
            Lsum = 0; move16();
            idx  = sub(c, k);
            FOR (; k < c; k++)
            {
                Lsum = L_add(Lsum, L_abs(x[ind[k]]));
            }
            fac_unq2 = BASOP_Util_Divide3216_Scale(Lsum, idx, &s1);
            fac_unq2 = BASOP_Util_Divide1616_Scale(fac_unq2, gg, &s2);
            s3       = sub(15, add(x_e, add(s1, sub(s2, gg_e))));
            s2       = norm_s(fac_unq1);
            test();
            IF (fac_unq2 != 0 && add(s3, s2) < 0)
            {
                fac_unq2 = MAX_16; move16();
            }
            ELSE
            {
                fac_unq2 = shr_r(fac_unq2, s_min(s3, 15));
            }

            /* calculate noise filling gain as minimum over high and low frequencies */
            fac_unq = s_min(fac_unq1, fac_unq2);
        }
        ELSE
        {
            Lsum = L_abs(x[ind[0]]);
            FOR (k = 1; k < c; k++)
            {
                Lsum = L_add(Lsum, L_abs(x[ind[k]]));
            }
            fac_unq = BASOP_Util_Divide3216_Scale(Lsum, c, &s1);
            fac_unq = BASOP_Util_Divide1616_Scale(fac_unq, gg, &s2);
            s3      = sub(15, add(x_e, add(s1, sub(s2, gg_e))));
            s2      = norm_s(fac_unq);
            test();
            IF (fac_unq != 0 && add(s3, s2) < 0)
            {
                fac_unq = MAX_16; move16();
            }
            ELSE
            {
                fac_unq = shr_r(fac_unq, s_min(s3, 15));
            }
        }
#else
        Lsum = L_abs(x[ind[0]]);
        FOR (k = 1; k < c; k++)
        {
            Lsum = L_add(Lsum, L_abs(x[ind[k]]));
        }
        fac_unq = BASOP_Util_Divide3216_Scale(Lsum, c, &s1);
        fac_unq = BASOP_Util_Divide1616_Scale(fac_unq, gg, &s2);
        s3      = sub(15, add(x_e, add(s1, sub(s2, gg_e))));
        s2      = norm_s(fac_unq);
        test();
        IF (fac_unq != 0 && add(s3, s2) < 0)
        {
            fac_unq = MAX_16; move16();
        }
        ELSE
        {
            fac_unq = shr_r(fac_unq, s_min(s3, 15));
        }

#endif
    }

    idx = round_fx(L_sub(0x80000, L_mult(fac_unq, 16)));
    if (sub(idx, 7) > 0)
    {
        idx = 7; move16();
    }
    if (idx < 0)
    {
        idx = 0; move16();
    }
    *fac_ns_idx = idx; move16();

    Dyn_Mem_Deluxe_Out();
}

