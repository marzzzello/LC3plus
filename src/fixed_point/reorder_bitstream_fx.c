/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"

#include "constants.h"
#include "functions.h"


void processReorderBitstream_fx(UWord8 *bytes, Word16 n_pccw, Word16 n_pc, Word16 b_left, Word8 *scratchBuffer)
{
    Word16        block_bytes;
    UWord8 *      bytes_tmp;

#ifdef DYNMEM_COUNT
    struct _dynmem
    {
        Word16        block_bits, block_bytes;
        UWord8 *      bytes_tmp;
    };
    Dyn_Mem_In("processReorderBitstream_fx", sizeof(struct _dynmem));
#endif

    bytes_tmp = (UWord8 *)scratchAlign(scratchBuffer, 0); /* Size = LC3_MAX_BYTES */

    if (n_pccw == 0)
    {
#ifdef DYNMEM_COUNT
        Dyn_Mem_Out();
#endif
        return;
    }

    assert(b_left >= 0);

    /* set block size in bits and full bytes */
    block_bytes = shr_sat(add(n_pc, 1), 1);

    /* rearrange bitstream */
    basop_memmove(&bytes_tmp[0], &bytes[b_left], block_bytes * sizeof(UWord8));
    basop_memmove(&bytes_tmp[block_bytes], &bytes[0], b_left * sizeof(UWord8));

    basop_memmove(&bytes[0], &bytes_tmp[0], add(block_bytes, b_left) * sizeof(UWord8));

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}


