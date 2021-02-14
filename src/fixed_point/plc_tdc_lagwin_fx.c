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


void processLagwin_fx(Word32 r[], const Word32 w[], Word16 m)
{
    /* Start Processing */
    Counter i;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processLagwin_fx", sizeof(struct { Counter i; }));
#endif

    FOR (i = 0; i < m; i++)
    {
        r[i + 1] = Mpy_32_32(r[i + 1], w[i]); move32();
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}


