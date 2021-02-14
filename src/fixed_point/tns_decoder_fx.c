/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "functions.h"


static Word32 IIRLattice(Word16 order, const Word16 *parCoeff, Word32 *state, Word32 x);

/*************************************************************************/

void processTnsDecoder_fx(Word16 rc_idx[], Word32 x[], Word16 xLen, Word16 order[], Word16 *x_e, Word16 BW_stopband_idx,
                          Word16 frame_dms, Word8 *scratchBuffer)
{
    Dyn_Mem_Deluxe_In(
        Word32 *state;
        Counter i, j;
        Word16  s1, s2, s, *rc, f, stopfreq, BW_stopband;
        Word16  numfilters, startfreq[TNS_NUMFILTERS_MAX];
    );

    state = (Word32 *)scratchAlign(scratchBuffer, 0);               /* Size = MAXLAG */
    rc    = (Word16 *)scratchAlign(state, sizeof(*state) * MAXLAG); /* Size = MAXLAG */

    numfilters  = 1;
    BW_stopband = BW_cutoff_bin_all[BW_stopband_idx]; move16();

    SWITCH (frame_dms)
    {
    case 25:
        startfreq[0] = 3; move16();
        BW_stopband  = shr_pos(BW_stopband, 2);
        BREAK;
    case 50:
        startfreq[0] = 6; move16();
        BW_stopband  = shr_pos(BW_stopband, 1);
        BREAK;
    case 100: startfreq[0] = 12; BREAK;
    }

    IF (sub(BW_stopband_idx, 3) >= 0 && frame_dms >= 50)
    {
        numfilters   = 2;
        startfreq[1] = shr_pos(BW_stopband, 1);
    }
    stopfreq = 0;

    test(); test();
    IF (order[0] > 0 || (sub(numfilters, 2) == 0 && order[1] > 0))
    {
        /* Scaling */
        f = startfreq[0]; move16();
        test();
        IF (sub(numfilters, 2) == 0 && order[0] == 0)
        {
            f = startfreq[1]; move16();
        }
        s1   = getScaleFactor32(x, f);
        s2   = getScaleFactor32(x + f, sub(xLen, f));
        s    = s_min(s1, sub(s2, 7)); /* 7 bits of headroom for IIR filtering */
        *x_e = sub(*x_e, s);

/* Init Filter */
        basop_memset(state, 0, MAXLAG * sizeof(Word32));
        FOR (i = 0; i < f; i++)
        {
            x[i] = L_shl(x[i], s); move32();
        }

        FOR (j = 0; j < numfilters; j++)
        {
            IF (order[j] > 0)
            {
                /* Unquantize coefficients */
                FOR (i = 0; i < order[j]; i++)
                {
                    rc[i] = tnsQuantPts[rc_idx[j * MAXLAG + i]]; move16();
                }

                /* Stop frequency */
                stopfreq = BW_stopband; move16();
                IF (sub(numfilters, 2) == 0 && j == 0)
                {
                    stopfreq = startfreq[1];
                }

                /* Filter */
                FOR (i = startfreq[j]; i < stopfreq; i++)
                {
                    x[i] = IIRLattice(order[j], rc, state, L_shl(x[i], s)); move32();
                }
            }
        }
        FOR (i = stopfreq; i < xLen; i++)
        {
            x[i] = L_shl(x[i], s); move32();
        }
    }
    Dyn_Mem_Deluxe_Out();
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

static Word32 IIRLattice(Word16 order, const Word16 *parCoeff, Word32 *state, Word32 x)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
    );

    /* first stage: no need to calculate state[order-1] */
    x = L_sub(x, Mpy_32_16(state[order - 1], parCoeff[order - 1]));

    FOR (i = order - 2; i >= 0; i--)
    {
        x            = L_sub(x, Mpy_32_16(state[i], parCoeff[i]));
        state[i + 1] = L_add(state[i], Mpy_32_16(x, parCoeff[i])); move32();
    }

    state[0] = x; move32();

    Dyn_Mem_Deluxe_Out();
    return x;
}

