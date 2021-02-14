/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"



void processResidualCoding_fx(Word16 x_e, Word32 x[], Word16 xq[], Word16 gain, Word16 gain_e, Word16 L_spec,
                              Word16 targetBits, Word16 nBits, Word8 *resBits, Word16 *numResBits)
{

    Counter i;
    Word16  s, n, m;
    Word32  L_tmp;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processResidualCoding_fx", sizeof(struct {
                   Counter i;
                   Word16  s, n, m;
                   Word32  L_tmp;
               }));
#endif

    n = 0; move16();
    m = add(sub(targetBits, nBits), 4);
    s = sub(add(15, gain_e), x_e);
    FOR (i = 0; i < L_spec; i++)
    {
        IF (xq[i] != 0)
        {
            L_tmp = L_sub(x[i], L_shl(L_mult(xq[i], gain), s));
            if (L_tmp < 0)
            {
                resBits[n] = 0; move16();
            }
            if (L_tmp >= 0)
            {
                resBits[n] = 1; move16();
            }
            n = add(n, 1);
            IF (sub(n, m) == 0)
            {
                BREAK;
            }
        }
    }
    *numResBits = n; move16();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

