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
st_PhECU_f0minLtp=55; % 55.4 Hz
st_PhECU_f0maxLtp=376; % 376.4706 Hz

*/


Word16
plc_phEcuSetF0Hz_fx(/*  output Q7  bin frequency [0.. 255.xxxx]  "1 sign, 8 bits mantissa, 7 binomial"  [0-255.9999]  */
                    Word16 fs_idx, Word16 old_pitch_int, Word16 old_pitch_fr)
{
    Word16 pitch_lagQ2, result, expo;
    Word32 L_result, L_tmp;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("plc_phEcuSetF0Hz_fx", sizeof(struct {
                   Word16 pitch_lagQ2, result, expo;
                   Word32 L_result, L_tmp;
                   Word16 num_FsByResQ0[5];
               }));
#endif

    BASOP_sub_sub_start("PhECU::plc_phEcuSetF0Hz_fx");

    result = 0; move16();
    IF (old_pitch_int != 0)
    {
        pitch_lagQ2 = add(
            old_pitch_fr,
            shl(old_pitch_int, 2)); /* lag at the current fs_idx , max lag_value is is 228(+.75)*48/12.8 = 858 in Q0 */

        L_result = plc_phEcu_ratio_fx(L_deposit_h(num_FsByResQ0[fs_idx]), L_deposit_h(pitch_lagQ2), &expo);
		L_tmp = L_shl_sat(L_result, sub(11, expo)); /* move to Q7, in high word to allow round*/
		result = round_fx(L_tmp);
    }
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
    BASOP_sub_sub_end();

    return result; /*Q7*/
}


