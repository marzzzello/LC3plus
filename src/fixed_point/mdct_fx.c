/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"



/* Union holding buffers to conserve stack memory. */

void processMdct_fx(
    Word16 x[],             /* i:   time input signal */
    Word16 x_exp, Word16 N, /* i:   block size N */
    const Word16 w[],       /* i:   window coefficients including normalization of sqrt(2/N) and scaled by 2^4 */
    Word16       wLen,      /* i:   window length */
    Word16       mem[],     /* i/o: last block of input samples */
    Word16       memLen,    /* i:   length of last sample block */
    Word32       y[],       /* o:   spectral data */
    Word16 *     y_e,       /* o:   spectal data exponent */
    Word8 *      scratchBuffer)
{
    Counter i;
    Word16  z, s, m;
    Word16 *buf;
    Word32 *workBuffer;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processMdct_fx", sizeof(struct {
                   Counter i;
                   Word16  z, s, m;
                   Word16 *buf;
                   Word32 *workBuffer;
               }));
#endif

    /* Buffers overlap since they are not used at the same time */
    buf        = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * MAX_LEN */
    workBuffer = (Word32 *)scratchAlign(scratchBuffer, 0); /* Size = 4 * MAX_LEN */

    /* Init (constant per sample rate) */
    z = (N << 1) - wLen; /* number of leading zeros in window */
    m = N >> 1;          /* half block size */

    basop_memmove(buf, mem, memLen * sizeof(Word16));

    basop_memmove(&buf[memLen], x, (N - memLen) * sizeof(Word16));

    basop_memmove(mem, &x[N - memLen], memLen * sizeof(Word16));

    FOR (i = 0; i < m; i++)
    {
        y[m + i] = L_msu0(L_mult0(buf[i], w[i]), buf[2 * m - 1 - i], w[2 * m - 1 - i]); move32();
    }

    FOR (i = 0; i < z; i++)
    {
        y[m - 1 - i] = L_mult0(x[2 * m - memLen + i], w[2 * m + i]); move32();
    }
    
    FOR (i = i; i < m; i++)
    {
        y[m - 1 - i] = L_mac0(L_mult0(x[2 * m - memLen + i], w[2 * m + i]), x[4 * m - memLen - 1 - i],
                              w[4 * m - 1 - i]); move32();
    }

    s = s_max(0, getScaleFactor32(y, N));
    FOR (i = 0; i < N; i++)
    {
        y[i] = L_shl(y[i], s); move32();
    }

    *y_e = sub(sub(x_exp, 2), s);

    /* N=20 only for 2.5ms possible */
    /* maybe implement this a pre init of shift */
    if (sub(N, 20) <= 0)
    {
        *y_e = add(*y_e, 2);
    }
    else if (sub(N, 120) <= 0)
    {
        *y_e = add(*y_e, 1);
    }

    dct_IV(y, y_e, N, workBuffer);

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

