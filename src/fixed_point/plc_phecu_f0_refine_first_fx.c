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


void plc_phEcu_F0_refine_first_fx(Word16 *     plocs,                       /* i/o */
                                  const Word16 n_plocs_in, Word32 *L_f0est, /* i/o  Q16 */
                                  const Word16 stPhECU_f0hzLtpBinQ7, const Word16 stPhECU_f0gainLtpQ15,
                                  const Word16 nSubm)

{
    Counter subm, i;
    Word16  ploc, n_plocs_ana;
    Word32  L_tmp = 0, L_diff, L_f0EstQ7, L_sensitivity_Q7;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("plc_phEcu_F0_refine_first_fx", sizeof(struct {
                   Counter subm, i;
                   Word16  ploc, n_plocs_ana;
                   Word32  L_tmp, L_diff, L_f0EstQ7, L_sensitivity_Q7;
               }));
#endif

     

    /* single initial peak F0 correction using available LTP information  */

    IF (sub(stPhECU_f0gainLtpQ15, ((Word16)(0.25 * 32768.0))) > 0)
    {
        ploc        = -1; move16();        /* sentinel */
        n_plocs_ana = s_min(n_plocs_in, 4); /* only analyze  at first 3 deteteced LF peaks */

        /*  only apply analysis below nsubm*pitmax_freq  ~=  1600Hz */
        i = sub(n_plocs_ana, 1);
        WHILE (i >= 0 && sub(plocs[i], (Word16)(1600.0 / 62.5)) > 0)
        {
            i--;
        }
        n_plocs_ana = add(i, 1);

        IF ((n_plocs_ana > 0))
        {
            /*   % find/correct first peak  in f0est , that is a submultiple of n*f0Ltp*/
            FOR (i = 0; i < n_plocs_ana; i++)
            {

                L_sensitivity_Q7 = L_deposit_l(((Word32)1) << (7 - 1)); /* 0.5 in Q7 */
                if (sub(stPhECU_f0gainLtpQ15, ((Word16)(0.75 * 32768.0))) < 0)
                {
                    L_sensitivity_Q7 = L_shr_pos(L_sensitivity_Q7, 1); /* %   more picky if correlation is rather low */
                }

                L_f0EstQ7 = L_shr_pos(L_f0est[i], 9); /* Q16 to Q7 */

                FOR (subm = 1; subm <= nSubm; subm++)
                {
                    /*adjf0 = abs(f0est - subm*stPhECU_f0hzLtpBin*ones(size(f0est))) < sensitivity ; % L1  difference,
                    vector operation over f0
                    ind   = find(adjf0==1,1); */
                    L_diff = L_msu0(L_f0EstQ7, subm, stPhECU_f0hzLtpBinQ7);
                    L_diff = L_abs(L_diff);
                    IF (L_sub(L_diff, L_sensitivity_Q7) < 0)
                    {
                        L_tmp = L_shl_pos(L_mult0(subm, stPhECU_f0hzLtpBinQ7), 16 - 7); /* to Q16 */
                        ploc  = i;                                                      move16();
                        BREAK;
                    }
                    L_sensitivity_Q7 = Mpy_32_16(L_sensitivity_Q7, (Word16)(0.875 * 32768.0 )); /* 2 cycles */
                }                                                                                    /* subm*/

                IF (ploc >= 0)
                {
                    BREAK;
                }
            } /* i, n_ploc_ana*/
        }

        if (ploc >= 0)
        {
            L_f0est[ploc] = L_tmp; move32(); /* in Q16 */
            /*ideally also integer plocs should be updated , e.g.  if  f0est goes from 1.45(plocs=1)  to 1.6(plocs==2)
             */
            /* but that is costly and   not required as long as corr_phase width is large enough ]*/
        }
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
   
}


