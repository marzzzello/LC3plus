/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"



void processApplyGlobalGain_fx(Word32 x[], Word16 *x_e, Word16 xLen, Word16 global_gain_idx, Word16 global_gain_off)
{
    Counter i;
    Word16 global_gain, global_gain_e;
    Word32 tmp32;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processApplyGlobalGain_fx", sizeof(struct {
                   Counter i;
                   Word16  global_gain, global_gain_e;
                   Word32  tmp32;
               }));
#endif

    tmp32         = L_shl_pos(L_mult0(add(global_gain_idx, global_gain_off), 0x797D), 7);
    global_gain_e = add(extract_l(L_shr_pos(tmp32, 25)), 1);
    global_gain = round_fx(BASOP_Util_InvLog2(L_or(tmp32, 0xFE000000)));

    FOR (i = 0; i < xLen; i++)
    {
        x[i] = Mpy_32_16(x[i], global_gain); move32();
    }

    *x_e = add(*x_e, global_gain_e); move16();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

