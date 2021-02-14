/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"
#include "functions.h"


/*
 * processTdac_fx
 *
 * Parameters:
 *   ola_mem       o: pointer of output signal                  Q0
 *   ola_mem_exp   o: exponent of output signal                 Q0
 *   synth         i: pointer of input signal                   Q0
 *   synth_exp     i: exponent of input signal                  Q0
 *   win           i: pointer of analysis and synthesis window  Q0
 *   la_zeroes     i: number of zeroes                          Q0
 *   frame_len     i: frame length                              Q0
 *
 * Function:
 *
 *
 * Returns:
 *    void
 */
void processTdac_fx(Word16 *ola_mem, Word16 *ola_mem_exp, const Word16 *synth_inp, const Word16 synth_exp_inp,
                    const Word16 *win, const Word16 la_zeroes, const Word16 frame_len, Word8 *scratchBuffer)
{
    Counter       i;
    Word16        s;
    Word16        L;
    Word16        N;
    Word16        NZ;
    Word16        LD2;
    Word32        sz;
    Word16        INV_NORM;
    Word16        INV_NORM_E;
    Word16        smax;
    Word16 *      synth;
    Word16        synth_len;
    Word16        synth_exp;
    const Word16 *win1;
    const Word16 *win2;
    const Word16 *win3;
    const Word16 *win4;
    const Word16 *synth1;
    const Word16 *synth2;
    Word16 *      ola_mem1;
    Word16 *      ola_mem2;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processTdac_fx", sizeof(struct {
                   Counter       i;
                   Word16        s;
                   Word16        L;
                   Word16        N;
                   Word16        NZ;
                   Word16        LD2;
                   Word32        sz;
                   Word16        INV_NORM;
                   Word16        INV_NORM_E;
                   Word16        smax;
                   Word16 *      synth;
                   Word16        synth_len;
                   Word16        synth_exp;
                   const Word16 *win1;
                   const Word16 *win2;
                   const Word16 *win3;
                   const Word16 *win4;
                   const Word16 *synth1;
                   const Word16 *synth2;
                   Word16 *      ola_mem1;
                   Word16 *      ola_mem2;
               }));
#endif

    synth = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * MAX_LEN */

    ASSERT(la_zeroes <= frame_len / 2);

    L   = frame_len; move16();
    LD2 = shr_pos(L, 1);
    NZ  = sub(LD2, la_zeroes);

    /* inverse normalization of sqrt(2/N) inside window */
    INV_NORM   = negate(shl_pos(frame_len, (15 - 9)));
    INV_NORM_E = 2; move16();
    if (norm_s(INV_NORM) > 0)
    {
        INV_NORM   = shl_pos(INV_NORM, 1);
        INV_NORM_E = 1; move16();
    }
    if (sub(frame_len, 120) <= 0)
    {
        INV_NORM_E = add(INV_NORM_E, 2);
    }
    if (sub(frame_len, 20) <= 0)
    {
        INV_NORM_E = add(INV_NORM_E, 2);
    }

    /* Scale input */
    synth_len = sub(shl_pos(L, 1), la_zeroes);
    s         = getScaleFactor16(synth_inp, synth_len);

    FOR (i = 0; i < synth_len; i++)
    {
        synth[i] = shl(synth_inp[i], s); move16();
    }
    synth_exp = sub(synth_exp_inp, s);

    /* calculate x_ov[L+la_zeroes] ... x_ov[2*L-1] */

    win1 = &win[L + LD2 - 1];
    win2 = &win[L + LD2];

    win3 = &win[LD2 - 1];
    win4 = &win[LD2];

    synth1 = &synth[L + LD2 - 1 - la_zeroes];
    synth2 = &synth[L + LD2 - la_zeroes];

    ola_mem1 = &ola_mem[LD2 - la_zeroes];
    ola_mem2 = &ola_mem[LD2 - la_zeroes - 1];

    smax = 15; move16();

    FOR (i = 0; i < NZ; i++)
    {
        /* analysis windowing + 2N -> N */
        sz = L_mac_sat(L_mult(*synth1, *win1), *synth2, *win2);

        /* N -> 2N + synthesis windowing */
        *ola_mem1 = round_fx(Mpy_32_16(sz, *win3)); move16();
        *ola_mem2 = round_fx(Mpy_32_16(sz, *win4)); move16();

        /* determine headroom */
        s = norm_s(*ola_mem1);
        if (*ola_mem1 != 0)
            smax = s_min(smax, s);
        s = norm_s(*ola_mem2);
        if (*ola_mem2 != 0)
            smax = s_min(smax, s);

        /* pointer update */
        win1--;
        win2++;
        win3--;
        win4++;
        synth1--;
        synth2++;
        ola_mem1++;
        ola_mem2--;
    }

    N = LD2; move16();

    FOR (; i < N; i++)
    {
        /* analysis windowing + 2N -> N */
        sz = L_mult(*synth1, *win1);

        /* N -> 2N + synthesis windowing */
        *ola_mem1 = round_fx(Mpy_32_16(sz, *win3)); move16();

        /* determin headroom */
        s = norm_s(*ola_mem1);
        if (*ola_mem1 != 0)
            smax = s_min(smax, s);

        /* pointer update */
        win1--;
        win2++;
        win3--;
        synth1--;
        synth2++;
        ola_mem1++;
    }

    smax = s_min(smax, 15);

    N = add(N, NZ);

    FOR (i = 0; i < N; i++)
    {
        ola_mem[i] = round_fx(L_mult(shl(ola_mem[i], smax), INV_NORM)); move16();
    }

    *ola_mem_exp = sub(add(synth_exp, INV_NORM_E), smax); move16();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}


