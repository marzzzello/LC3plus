/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/


#include "functions.h"


void dct_IV(Word32 *pDat,       /* i/o: pointer to data buffer */
            Word16 *pDat_e,     /* i/o: pointer to data exponent */
            Word16  L,          /* i  : length of block */
            Word32 *workBuffer) /* : size of L */
{
    Word16 sin_step;
    Word16 idx;
    Word16 M_var;
    Word16 M2;

    Word32 *pDat_0;
    Word32 *pDat_1;

    Word32 accu1;
    Word32 accu2;
    Word32 accu3;
    Word32 accu4;

    Counter i;

    const PWord16 *twiddle;
    const PWord16 *sin_twiddle;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("dct_IV", sizeof(struct {
                   Word16  sin_step;
                   Word16  idx;
                   Counter i;
                   Word16  M_var;
                   Word16  M2;

                   Word32 *pDat_0;
                   Word32 *pDat_1;

                   Word32 accu1;
                   Word32 accu2;
                   Word32 accu3;
                   Word32 accu4;

                   const PWord16 *twiddle;
                   const PWord16 *sin_twiddle;
               }));
#endif

    M_var = shr_pos_pos(L, 1);
    M2    = sub(M_var, 1);

    BASOP_getTables(&twiddle, &sin_twiddle, &sin_step, L);

    pDat_0 = &pDat[0];
    pDat_1 = &pDat[L - 2];

    FOR (i = 0; i < M2; i += 2)
    {
        cplxMpy32_32_16_2(accu1, accu2, pDat_1[1], pDat_0[0], twiddle[i].v.re, twiddle[i].v.im);
        cplxMpy32_32_16_2(accu3, accu4, pDat_1[0], pDat_0[1], twiddle[i + 1].v.re, twiddle[i + 1].v.im);

        pDat_0[0] = accu2;           move32();
        pDat_0[1] = accu1;           move32();
        pDat_1[0] = accu4;           move32();
        pDat_1[1] = L_negate(accu3); move32();

        pDat_0 = pDat_0 + 2;
        pDat_1 = pDat_1 - 2;
    }

    BASOP_cfft(&pDat[0], &pDat[1], M_var, 2, pDat_e, workBuffer);

    pDat_0 = &pDat[0];
    pDat_1 = &pDat[L - 2];

    idx = sin_step;
    M2  = sub(shr_pos_pos(add(M_var, 1), 1), 1);

    /* Sin and Cos values are 0.0f and 1.0f */
    cplxMpy32_32_16_2(accu3, accu4, pDat_1[0], pDat_1[1], sin_twiddle[idx].v.re, sin_twiddle[idx].v.im);

    pDat_1[1] = L_negate(L_shr_pos(pDat_0[1], 1)); move32();
    pDat_0[0] = L_shr_pos(pDat_0[0], 1);           move32();

    FOR (i = 1; i < M2; i++)
    {
        pDat_0[1] = accu3; move32();
        pDat_1[0] = accu4; move32();

        pDat_0 = pDat_0 + 2;
        pDat_1 = pDat_1 - 2;

        cplxMpy32_32_16_2(accu1, accu2, pDat_0[1], pDat_0[0], sin_twiddle[idx].v.re, sin_twiddle[idx].v.im);

        idx += sin_step;

        cplxMpy32_32_16_2(accu3, accu4, pDat_1[0], pDat_1[1], sin_twiddle[idx].v.re, sin_twiddle[idx].v.im);

        pDat_1[1] = L_negate(accu1); move32();
        pDat_0[0] = accu2;           move32();
    }

    pDat_0[1] = accu3; move32();
    pDat_1[0] = accu4; move32();

    pDat_0 = pDat_0 + 2;
    pDat_1 = pDat_1 - 2;

    cplxMpy32_32_16_2(accu3, accu4, pDat_0[1], pDat_0[0], sin_twiddle[idx].v.re, sin_twiddle[idx].v.im);

/* Last Sin and Cos value pair are the same */
    accu1  = L_shr_pos(Mpy_32_16(pDat_1[0], TWIDDLE), 1);
    accu2  = L_shr_pos(Mpy_32_16(pDat_1[1], TWIDDLE), 1);

    pDat_1[0] = L_add(accu1, accu2); move32();
    pDat_0[1] = L_sub(accu1, accu2); move32();

    pDat_1[1] = L_negate(accu3); move32();
    pDat_0[0] = accu4;           move32();

    /* twiddeling scale is 2 */
    *pDat_e = add(*pDat_e, 2); move16();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

