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


void plc_phEcu_LF_peak_analysis_fx(Word16 *      plocs,      /* i/o  0 ... Lprot/2 +1*/
                                   Word16 *      n_plocs,    /* i/o   0.. MAX_PLOCS  */
                                   Word32 *      L_f0estQ16, /* i/o  Q16*/
                                   const Word16 *mag,        /* i: Qx    */
                                   const Word16 stPhECU_f0hzLtpBinQ7, const Word16 stPhECU_f0gainLtpQ15,
                                   const Word16 nSubm, Word16 maxPlocs,
                                   Word8 *scratchBuffer /* Size = 6 * MAX_PLOCS + 42 */
)

{
    Counter i, j;
    Word16  n_plocs_ana, peakLF_Xval, tmp, f_abs_ind, plocsIntersectFlag;

    Word32  L_fQ7, *L_f0est_prelQ16;
    Word16  num_prel = 0, *plocs_prel;
    Word16  prel_low, prel_high, start, fin;
    Word16 *plocs_old;
    Word32 *L_plocsi_old;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("plc_phEcu_LF_peak_analysis_fx", sizeof(struct {
                   Counter i, j;
                   Word16  n_plocs_ana, peakLF_Xval, tmp, f_abs_ind, plocsIntersectFlag;
                   Word32  L_fQ7, *L_f0est_prelQ16;
                   Word16  num_prel, *plocs_prel;
                   Word16  prel_low, prel_high, start, fin;
                   Word16 *plocs_old;
                   Word32 *L_plocsi_old;
               }));
#endif

     

    L_f0est_prelQ16 = (Word32 *)scratchAlign(scratchBuffer, 0);                              /* Size = 4 * 7 */
    plocs_prel      = (Word16 *)scratchAlign(L_f0est_prelQ16, sizeof(*L_f0est_prelQ16) * 7); /* Size = 2 * 7 */
    plocs_old       = (Word16 *)scratchAlign(plocs_prel, sizeof(*plocs_prel) * 7);           /* Size = 2 * MAX_PLOCS */
    L_plocsi_old    = (Word32 *)scratchAlign(plocs_old, sizeof(*plocs_old) * maxPlocs);      /* Size = 4 * MAX_PLOCS */

    test(); test();
    IF ((*n_plocs > 0) && sub(stPhECU_f0gainLtpQ15, ((Word16)(0.25 * 32768.0))) > 0 &&
        sub(stPhECU_f0hzLtpBinQ7, (Word16)(2.75 * 128.0)) < 0)
    {

        /* % analyze/apply  f0Ltp to avoid  intermodulation effects  below  f0  of ~180 Hz
        % we only do the  f0Ltp-replacement(s)  if  there is already an established
        % fft peak in the region   ~fRes  to  2.5*fres
        fft_peak_eval_plocs = 1:3;
        plocsIntersectFlag = intersect(plocs, fft_peak_eval_plocs );  % check for 1,2,3  in plocs  */

        plocsIntersectFlag = 0; move16();
        peakLF_Xval        = 0; move16();
        n_plocs_ana        = s_min(*n_plocs, 3);
        FOR (i = 0; i < n_plocs_ana; i++)
        {
            tmp = plocs[i];       move16();
            if (sub(tmp, 2) <= 0) /*  C index  0, 1,2  checked , [DC, 62.5 Hz, 125Hz ] */
            {
                plocsIntersectFlag = add(i, 1);
            }
            peakLF_Xval = s_max(mag[tmp], peakLF_Xval);
        }

        num_prel = 0; move16();
        IF (plocsIntersectFlag != 0)
        { /* fft-peak at 0, 62 or 125 Hz  */
            /*  analyze if  ltp-based f0 need to be added  or not  */
            peakLF_Xval = mult_r(peakLF_Xval, (Word16)(.375 * 32768.0)); /* now as a limit */

            FOR (i = 1; i <= nSubm; i++)
            {
                L_fQ7     = L_mult0(i, stPhECU_f0hzLtpBinQ7); /* fractional index stored in L_plocsi */
                f_abs_ind = L_shr_pos(L_add(L_fQ7, 64), 7);   /* integer bin index stored in plocs */

                test();
                IF ((L_sub(L_fQ7, 819) <= 0) && /*  % only apply up to ~400hz , 819 = 400/62.5*128 */
                    (sub(mag[f_abs_ind], peakLF_Xval) >
                     0)) /* %  only set as preliminary  if relative peak strength is signficant*/
                {
                    L_f0est_prelQ16[num_prel] = L_shl_pos(L_fQ7, 9); move32();
                    plocs_prel[num_prel]      = f_abs_ind;           move16();
                    num_prel                  = add(num_prel, 1);
                }
            }
        } /*intersectFlag*/

        /* now replace/ merge new preliminary added peaks with existing plocs and L_f0estQ16 */
        /* note that a previous fake/merged magnitude-determined peak may be replaced by two separated  side peaks */

        /* a general non-optimized list-merging solution below */
        test();
        IF ((num_prel > 0) && (sub(add(num_prel, *n_plocs), MAX_PLOCS) <= 0) /* skip in case plocs list is too large */
        )
        {
            prel_low  = plocs_prel[0];                move16();
            prel_high = plocs_prel[sub(num_prel, 1)]; move16();

            start = -1; move16();
            FOR (i = sub(*n_plocs, 1); i >= 0; i--)
            {
                if (sub(plocs[i], prel_low) >= 0)
                {
                    start = i; move16();
                }
            }
            start = sub(start, 1);    /* end of old section to copy before the added/merged section */
            start = s_max(start, -1); /* limit  for loop later */
                                      /*% dbg check  low part for a sucessful replace/merge  */
            if (start >= 0 && start < *n_plocs)
            {
                ASSERT(plocs[start] < plocs_prel[0]);
            }

            sub(0, 0);
            IF (prel_high < plocs[0])
            {
                fin = 0; move16(); /*% keep all plocs , just concat  */
            }
            ELSE
            {
                fin = *n_plocs;
                FOR (i = 0; i < *n_plocs; i++)
                {
                    sub(0, 0);
                    if (plocs[i] <= prel_high)
                    {
                        fin = i; move16();
                    }
                }
                fin = add(fin, 1); /* first element  in high part of old  plocs to be copied  */
            }

            /*% dbg check high part for a sucessful replace/merge */
            if (fin >= 0 && fin < *n_plocs)
            {
                ASSERT(plocs_prel[sub(num_prel, 1)] < plocs[fin]);
            }

            /*
            % actual replace/merge of added integer locations and fractional freqs. into plocs/f0list  list ;
            % three loops in BASOP
            plocs     =  [ plocs(1:(start)) ; plocs_prel ; plocs((fin):end) ];
            f0est    =   [  f0est(1:(start)) ; f0est_prel; f0est((fin):end) ];
            */

            FOR (i = 0; i < *n_plocs; i++)
            {
                plocs_old[i]    = plocs[i];      move16();
                L_plocsi_old[i] = L_f0estQ16[i]; move32();
            }

            /*
            j=0;
            FOR(i=0; i <= start; i++)
            {
                plocs[i] = plocs_old[i];        move16();
            L_f0estQ16[i] = L_plocsi_old[i]; move32();
            j++;
        }
        */

            j = add(start, 1);

        FOR (i = 0; i < num_prel; i++) /* NB this section may  both insert or overwrite old plocs   */
        {
            plocs[j]      = plocs_prel[i];      move16();
            L_f0estQ16[j] = L_f0est_prelQ16[i]; move32();
            j++;
        }
        FOR (i = fin; i < *n_plocs; i++) /* copy the tail of the list */
        {
            plocs[j]      = plocs_old[i];    move16();
            L_f0estQ16[j] = L_plocsi_old[i]; move32();
            j++;
        }

        *n_plocs = j; move16(); /* update total length   */
    }                            /* num_prel >0*/
} /* gain/hz Limits */

#ifdef DYNMEM_COUNT
Dyn_Mem_Out();
#endif
    
}


