/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"



void processResidualDecoding_fx(Word32 x[], Word16 x_e, Word16 L_spec, Word16 prm[], Word16 resQBits)
{

    Counter i;
    Word32  fac_m, fac_p;
    Word16  s, bits;
    Word32  tmp;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processResidualDecoding_fx", sizeof(struct {
                   Counter i;
                   Word32  fac_m, fac_p;
                   Word16  s, bits;
                   Word32  tmp;
               }));
#endif

    tmp   = 0;
    s     = sub(x_e, 1);
    fac_m = L_shr(0xC000000, s);
    fac_p = L_shr(0x14000000, s);
    bits  = 0; move16();

    FOR (i = 0; i < L_spec; i++)
    {
        IF (sub(bits, resQBits) >= 0)
        {
            BREAK;
        }

        IF (x[i] != 0)
        {
            IF (prm[bits] == 0)
            {
                if (x[i] > 0)
                    tmp = L_sub(x[i], fac_m);
                if (x[i] < 0)
                    tmp = L_sub(x[i], fac_p);
            }
            ELSE
            {
                if (x[i] > 0)
                    tmp = L_add(x[i], fac_p);
                if (x[i] < 0)
                    tmp = L_add(x[i], fac_m);
            }
            x[i] = tmp; move32();
            bits = add(bits, 1);
        }
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

