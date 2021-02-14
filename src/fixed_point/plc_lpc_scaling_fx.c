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


void processPLCLpcScaling_fx(Word32 tdc_A_32[], Word16 tdc_A_16[], Word16 m)
{
    Counter i;
    Word16  s;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processPLCLpcScaling_fx", sizeof(struct {
                   Counter i;
                   Word16  s;
               }));
#endif

    s = getScaleFactor32(tdc_A_32, m);
    FOR (i = 0; i < m; i++)
    {
        tdc_A_16[i] = round_fx_sat(L_shl_sat(tdc_A_32[i], s)); move16();
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}


