/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "basop_mpy.h"
#include "stl.h"


Word32 Mpy_32_16(Word32 x, Word16 y)
{
    Word32  mh;
    UWord16 ml;

    Mpy_32_16_ss(x, y, &mh, &ml);

    return (mh);
}


Word32 Mpy_32_32(Word32 x, Word32 y)
{
    Word32  mh;
    UWord32 ml;

    Mpy_32_32_ss(x, y, &mh, &ml);

    return (mh);
}


void cplxMpy_32_16(Word32 *c_Re, Word32 *c_Im, const Word32 a_Re, const Word32 a_Im, const Word16 b_Re,
                   const Word16 b_Im)
{
    *c_Re = L_sub(Mpy_32_16(a_Re, b_Re), Mpy_32_16(a_Im, b_Im)); move32();
    *c_Im = L_add(Mpy_32_16(a_Re, b_Im), Mpy_32_16(a_Im, b_Re)); move32();
}

