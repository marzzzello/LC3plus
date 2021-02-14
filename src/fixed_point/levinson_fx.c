/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"



void processLevinson_fx(Word32 *lpc, Word32 *ac, Word16 N, Word16 *rc, Word32 *pred_err, Word8 *scratchBuffer)
{

    Word32 *lpc_tmp;
    Word32  rc32, err, sum;
    Word16  shift, s, inv;
    Counter n, m;


#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processLevinson_fx", sizeof(struct {
                   Word32 *lpc_tmp;
                   Word32  rc32, err, sum;
                   Word16  shift, s, inv;
                   Counter n, m;
                   Word32  params[2];
               }));
#endif

    lpc_tmp = (Word32 *)scratchAlign(scratchBuffer, 0); /* Size = 4 * (M_LTPF + 1) = 100 bytes */

    /* Init Prediction Error */
    err   = ac[0]; move32();
    shift = 0;     move16();

    /* LPC Coefficient 0 */
    lpc[0] = 0x8000000; move32();

    /* Reflection Coefficient 0 */
    IF (ac[0] != 0)
    {
        inv  = div_s(16383, extract_h(ac[0]));
        rc32 = L_shl_pos(Mpy_32_32(L_abs(ac[1]), Mpy_32_16(L_sub(MAX_32, Mpy_32_16(ac[0], inv)), inv)), 2);
    }
    ELSE
    {
        rc32 = 0; move32();
    }
    if (ac[1] > 0)
    {
        rc32 = L_negate(rc32);
    }
    if (rc != NULL)
    {
        rc[0] = round_fx(rc32); move16();
    }

    /* LPC Coefficient 1 */
    lpc[1] = L_shr_pos(rc32, 4); move32();

    FOR (n = 2; n <= N; n++)
    {
        /* Update Prediction Error */
        err   = Mpy_32_32(err, L_sub(MAX_32, Mpy_32_32(rc32, rc32)));
        s     = norm_l(err);
        err   = L_shl_pos(err, s);
        shift = add(shift, s);

        /* Reflection Coefficient n-1 */
        sum = Mpy_32_32(ac[1], lpc[n - 1]);
        FOR (m = 2; m < n; m++)
        {
            sum = L_add(sum, Mpy_32_32(ac[m], lpc[n - m]));
        }

        sum = L_add(L_shl_pos(sum, 4), ac[n]);
        IF (err != 0)
        {
            inv  = div_s(16383, extract_h(err));
            rc32 = L_shl_pos(Mpy_32_32(L_abs(sum), Mpy_32_16(L_sub(MAX_32, Mpy_32_16(err, inv)), inv)), 2);
        }
        ELSE
        {
            rc32 = 0;
        }
        if (sum > 0)
        {
            rc32 = L_negate(rc32);
        }
        rc32 = L_shl(rc32, shift);
        if (rc != NULL)
        {
            rc[n - 1] = round_fx(rc32); move16();
        }

/* Recompute LPC Coefficients up to n-1 */
        FOR (m = 1; m < n; m++)
        {
            lpc_tmp[m] = L_add(lpc[m], Mpy_32_32(rc32, lpc[n - m])); move32();
        }

        basop_memmove(&lpc[1], &lpc_tmp[1], (n - 1) * sizeof(Word32));

        /* LPC Coefficient n */
        lpc[n] = L_shr_pos(rc32, 4); move32();
    }

    /* Final Prediction Error */
    IF (pred_err != NULL)
    {
        err       = Mpy_32_32(err, L_sub(MAX_32, Mpy_32_32(rc32, rc32)));
        *pred_err = L_shr(err, shift);
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}


void lpc2rc(Word32 *lpc, Word16 *rc, Word16 N)
{
    Word32  lpc_tmp[MAXLAG + 1];
    Word32  rc32, tmp0, tmp1;
    Word16  inv;
    Counter n, m;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("lpc2rc", sizeof(struct {
                   Word32  lpc_tmp[MAXLAG + 1];
                   Word32  rc32, tmp0, tmp1;
                   Word16  inv;
                   Counter n, m;
               }));
#endif

    FOR (n = N; n >= 2; n--)
    {
        rc32      = L_shl_pos(lpc[n], 4);
        rc[n - 1] = round_fx(rc32); move16();

        tmp0 = L_sub(MAX_32, L_abs(Mpy_32_32(rc32, rc32)));
        FOR (m = 1; m < n; m++)
        {
            tmp1       = L_sub(lpc[m], Mpy_32_32(lpc[n - m], rc32));
            inv        = div_s(16383, extract_h(tmp0));
            lpc_tmp[m] = L_shl_pos(Mpy_32_32(tmp1, Mpy_32_16(L_sub(MAX_32, Mpy_32_16(tmp0, inv)), inv)), 2); move32();
        }

        basop_memmove(&lpc[1], &lpc_tmp[1], (n - 1) * sizeof(Word32));
    }

    rc[0] = round_fx(L_shl_pos(lpc[1], 4)); move32();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

