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




void processPLCUpdateAfterIMDCT_fx(Word16 x_fx[], Word16 q_fx_exp, Word16 concealMethod, Word16 xLen, Word16 fs_idx,
   Word16 *nbLostFramesInRow, Word16 *prev_prev_bfi, Word16 *prev_bfi, Word16 bfi, Word16 scf_q[],
   Word16 *ns_cum_alpha, AplcSetup *plcAd)
{
#ifdef BE_MOVED_STAB_FAC
   Word16  oldLen, bufHistlen;
   Word16  scale_fac_old, scale_fac_new, q_theo_new_old, q_theo_new_new, q_new, shift_old, shift_new;
   Word16 frontLen, pastLen;
   Word16 marginOldPast;
   Word16 scale_fac_old_dual;
   Word16   marginNewXlen, marginOldFront;
#else
   Word16  oldLen, bufHistlen, d;
   Counter i;
   Word16  scale_fac_old, scale_fac_new, q_theo_new_old, q_theo_new_new, q_new, shift_old, shift_new;
   Word16 frontLen, pastLen;
   Word16 marginOldPast;
   Word16 scale_fac_old_dual;
   Word32  tmp32;
   Word16   marginNewXlen, marginOldFront;
#endif

#ifdef DYNMEM_COUNT
   Dyn_Mem_In("processPLCUpdateAfterIMDCT_fx", sizeof(struct {
      Word16  oldLen, bufHistlen;
      Counter i;
      Word16 scale_fac_old, scale_fac_new, q_theo_new_old, q_theo_new_new, q_new, shift_old, shift_new;
      Word16 frontLen, pastLen;
      Word16 marginOldPast;
      Word16 scale_fac_old_dual;
      Word16   marginNewXlen, marginOldFront;
   }));
#endif



    
      BASOP_sub_sub_start("processPLCUpdateAfterIMDCT ");
   
  
   if (plcAd)
   {
#ifdef       NONBE_FIX_PCMHIST_LENGTHS
      /* for  short NB frames(2.5 ms)  TDC-filtering  requires  more PCM samples than  the plc_xcorr function */
      bufHistlen = s_max(xLen,  add(M + 1, shr(xLen, 1))) ; 

      bufHistlen = add(pitch_max[fs_idx], bufHistlen );       
#else
      bufHistlen = add(pitch_max[fs_idx], xLen);
#endif


      logic16();
      IF( (sub(bfi,1)== 0 )  && sub(concealMethod, 2) == 0)
      {   /* % reduced buffering update length during concealment method 2 as Xsav_fx is stored in the  joint  q_old_fx and x_old_tot_fx buffer */
         bufHistlen = sub(bufHistlen, sub(LprotSzPtr[fs_idx], s_min(MAX_BW_BIN, xLen)));
         ASSERT(xLen == (Word16)(((double)LprotSzPtr[fs_idx])*0.625)); /*/ only enter here for 10 ms cases */

         /* actually one can  select to always update xLen(10 ms)  less  samples of x_old_tot,  also in  TDC-PLC bfi frames ,, and for PhECU.PLC  */
      }
      oldLen = sub(bufHistlen, xLen);

      /* update ltpf-free pcm history buffer for TD-PLC */

      basop_memmove(&plcAd->x_old_tot_fx[plcAd->max_len_pcm_plc - bufHistlen],
         &plcAd->x_old_tot_fx[plcAd->max_len_pcm_plc - bufHistlen + xLen], oldLen * sizeof(Word16));

      basop_memcpy(&plcAd->x_old_tot_fx[plcAd->max_len_pcm_plc - xLen], &x_fx[0], xLen * sizeof(Word16));

      frontLen = sub(LprotSzPtr[fs_idx], xLen);  /*16-10 =  6ms  of the  prev_synth/xfp part  */
      pastLen = sub(oldLen, frontLen);          /* ~11.8 ms*/

      marginOldPast = getScaleFactor16_0(&(plcAd->x_old_tot_fx[plcAd->max_len_pcm_plc - bufHistlen]), pastLen);
      marginOldFront = getScaleFactor16_0(&(plcAd->x_old_tot_fx[plcAd->max_len_pcm_plc - bufHistlen + pastLen]), frontLen);

      scale_fac_old_dual = s_min(marginOldFront, marginOldPast);
      scale_fac_old = scale_fac_old_dual;
      
      frontLen = 0; move16();
      logic16();  logic16();
      IF(bfi == 1 && *prev_bfi == 0 && sub(concealMethod, 2) == 0)
      {   /* prepare localized margin_xfp value  for a next bad concealment Method 2 frame   */
         frontLen = *nbLostFramesInRow;
         frontLen = add(hamm_len2Tab[fs_idx], shr(hamm_len2Tab[fs_idx], 2)); /*  find margin in the   3.75 ms front part   */
         pastLen = sub(xLen, frontLen);
         scale_fac_new = getScaleFactor16_0(&(x_fx[0]), pastLen);
         marginNewXlen = getScaleFactor16_0(&(x_fx[0]) + pastLen, frontLen); /* for pHEcuprev_synth  in 2nd+  bfi frame */

         scale_fac_new = s_min(scale_fac_new, marginNewXlen);
      }
      ELSE
      { /* prepare margin value for any coming  good frame  or  any coming first bad frame  */

          marginNewXlen = getScaleFactor16_0(&(x_fx[0]),xLen);  /* prevsynth  in first bfi frame */
          scale_fac_new = marginNewXlen; move16();
      }

      q_theo_new_old = s_max(plcAd->q_fx_old_exp - scale_fac_old, 0);
      q_theo_new_new = s_max(q_fx_exp - scale_fac_new, 0);

      q_new = s_max(q_theo_new_old, q_theo_new_new);

      shift_old = plcAd->q_fx_old_exp - q_new;
      shift_new = q_fx_exp - q_new;

      IF(shift_old != 0)
      {
         Scale_sig(&plcAd->x_old_tot_fx[plcAd->max_len_pcm_plc - bufHistlen], oldLen, shift_old);
         logic16();
         if ((bfi == 1) && (sub(concealMethod, 3) == 0))
         {
            plcAd->tdc_gain_c = L_shl(plcAd->tdc_gain_c, shift_old);
         }
         move16(); /* count move to static RAM */

         marginOldFront = s_min(16, sub(marginOldFront, shift_old));
      }
      IF(shift_new)
      {
         Scale_sig(&plcAd->x_old_tot_fx[plcAd->max_len_pcm_plc - xLen], xLen, shift_new); /* positive shift_new means upshift=less margin  */
         marginNewXlen = s_min(16, sub(marginNewXlen, shift_new));
      }

      plcAd->q_fx_old_exp = sub(q_fx_exp, shift_new);

      plcAd->PhECU_margin_xfp = s_min(marginNewXlen, marginOldFront);  move16(); /* for pHECU winEncalc xfp energy calculations */
      if (frontLen != 0)
      {  /* prepare margin value for a first pHECU(16 ms)  or a consecutive bad PhEcu frame (3.75ms)  */
         plcAd->PhECU_margin_xfp = marginNewXlen; move16();
      }
      if (sub(plcAd->PhECU_margin_xfp, 16) == 0)
      {
         plcAd->PhECU_margin_xfp = 1;   move16();  /* "1" --> does not rescale the   all-zero vector, inside PhECU  */
      }
   }

   /* Update PLC params */
   IF(sub(bfi, 1) != 0)
   {
      /* % reset counters in GF  */
      *nbLostFramesInRow = 0;  move16();  /*plc0,3 4 , udpate  */
      *ns_cum_alpha = 32767;  move16();   /*plc0,  4 , udpate  */

      if (plcAd)
      {
#ifndef BE_MOVED_STAB_FAC
         /* calculate stability factor */
         IF(sub(*prev_bfi, 1) == 0)
         {
            plcAd->stab_fac = 26214;  move16();
         }
         ELSE
         {
             tmp32 = 0;  move32();
             FOR(i = 0; i < M; i++)
             {
                 d = sub(scf_q[i], plcAd->old_scf_q[i]);
                 tmp32 = L_mac_sat(tmp32, d, d);
             }
             tmp32 = L_shl_sat(tmp32, 3);
             IF(tmp32 > 0x7D000000 /*1.25*25*/)
             {
                 plcAd->stab_fac = 0;  move16();
             }
             ELSE IF(tmp32 < 0x19003E82 /*0.25*25*/)
             {
                 plcAd->stab_fac = 0x7FFF;  move16();
             }
             ELSE
             {
                 tmp32 = L_shl_pos(L_sub(0x50000000 /*1.25/2*/, Mpy_32_16(tmp32, 0x51EC /*16/25*/)), 1);
                 plcAd->stab_fac = round_fx(tmp32);  move16();
             }
         }
#endif

         basop_memmove(plcAd->old_old_scf_q, plcAd->old_scf_q, M * sizeof(Word16));
         basop_memmove(plcAd->old_scf_q, scf_q, M * sizeof(Word16));

         /* PLC fullband transient detector setting for non-bfi frames */
         plcAd->PhECU_short_flag_prev = 0;  move16(); /* fullband transient not active   */
      }
   }
 
   /* values may be {0,1,2} */
   *prev_prev_bfi = *prev_bfi;  move16();
   *prev_bfi = bfi;  move16();

    BASOP_sub_sub_end();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

#ifdef BE_MOVED_STAB_FAC
void processPLCcomputeStabFac_main(Word16 scf_q[], Word16 old_scf_q[], Word16 old_old_scf_q[], Word16 bfi, Word16 prev_bfi,
                              Word16 prev_prev_bfi, Word16 *stab_fac)
{
    IF (sub(bfi, 1) == 0)
    {
        IF (sub(prev_bfi, 1) != 0)
        {
            processPLCcomputeStabFac(old_scf_q, old_old_scf_q, prev_prev_bfi, stab_fac);
        }
    }
    ELSE IF(sub(bfi, 2) == 0)
    {
        processPLCcomputeStabFac(scf_q, old_scf_q, prev_bfi, stab_fac);
    }
}

void processPLCcomputeStabFac(Word16 scf_q[], Word16 old_scf_q[], Word16 prev_bfi, Word16 *stab_fac)
{
    Counter i;
    Word32  tmp32;
    Word16  d;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("calculateStabFac", sizeof(struct {
                   Counter i;
                   Word32  tmp32;
                   Word16  d;
               }));
#endif

    /* calculate stability factor */
    IF (sub(prev_bfi, 1) == 0)
    {
        *stab_fac = 26214; move16();
    }
    ELSE
    {
        tmp32 = 0; move32();
        FOR (i = 0; i < M; i++)
        {
            d     = sub(scf_q[i], old_scf_q[i]);
            tmp32 = L_mac_sat(tmp32, d, d);
        }
        tmp32 = L_shl_sat(tmp32, 3);
        IF (tmp32 > 0x7D000000 /*1.25*25*/)
        {
            *stab_fac = 0; move16();
        }
        ELSE IF (tmp32 < 0x19003E82 /*0.25*25*/)
        {
            *stab_fac = 0x7FFF; move16();
        }
        ELSE
        {
            tmp32     = L_shl_pos(L_sub(0x50000000 /*1.25/2*/, Mpy_32_16(tmp32, 0x51EC /*16/25*/)), 1);
            *stab_fac = round_fx(tmp32); move16();
        }
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}
#endif /* BE_MOVED_STAB_FAC */

void processPLCUpdateXFP_w_E_hist_fx(Word16 prev_bfi, Word16 bfi, Word16 *xfp_fx, Word16 xfp_exp_fx, Word16 margin_xfp, 
                                     Word16 fs_idx,
                                     Word32 *L_oold_xfp_w_E_fx, Word16 *oold_xfp_w_E_exp_fx, 
                                     Word32 *L_old_xfp_w_E_fx, Word16 *old_xfp_w_E_exp_fx,
                                    
                                     Word16 *oold_Ltot_exp_fx ,Word16 *old_Ltot_exp_fx )    

{
    Word32 L_tot  ; 
    Word16 dn_scale, exp_shift;
    Word16 used_xfp_exp_fx;
    Word16 exp_out  ; 

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("PLCUpdateXFP_w_E_hist", sizeof(struct {  
        Word32 L_tot;
        Word16 dn_scale, exp_shift;
        Word16 used_xfp_exp_fx;
        Word16 exp_out; 
    }));
#endif
    BASOP_sub_sub_start("PhECU::UpdateXfp_w_E_hist_fx");

    IF (sub(bfi,1) != 0)
    {

        if (sub(prev_bfi,1) == 0)
        {
          /* only a single historic frame available in the next  frame  
               , force artifical update of oold energy to be the same as old */
           *old_xfp_w_E_exp_fx = LTOT_INIT_FLAG ;  move16();
        }

        /* Time shift energy state and xfp exp */ 
        IF (sub_sat(*old_xfp_w_E_exp_fx, LTOT_INIT_FLAG ) ==  0) 
        {
            *L_oold_xfp_w_E_fx   =   LTOT_MIN_MAN  ;                             move32();   
            *oold_xfp_w_E_exp_fx =  UNINIT_OR_UNSAFE_OOLD_SENTINEL ; move16();  
        }
        ELSE
        {
            *L_oold_xfp_w_E_fx   = *L_old_xfp_w_E_fx;   move32();  /* regular update */
            *oold_xfp_w_E_exp_fx = *old_xfp_w_E_exp_fx; move16();
        }

        /* Time shift L_tot energy state and L_tot_exp  */
        IF (sub_sat(*old_Ltot_exp_fx, LTOT_INIT_FLAG ) ==  0) 
        {
            *L_oold_xfp_w_E_fx   =   LTOT_MIN_MAN  ;                             move32();   
            *oold_Ltot_exp_fx    =   UNINIT_OR_UNSAFE_OOLD_SENTINEL ;    move16();  
        }
        ELSE
        {
            *L_oold_xfp_w_E_fx   = *L_old_xfp_w_E_fx;   move32();  /* regular update */
            *oold_Ltot_exp_fx    = *old_Ltot_exp_fx;        move16();
        }

       
        dn_scale        = e_tot_headroom[fs_idx]; /* allowed minimum dn_scale for a max upshifted signal */
        used_xfp_exp_fx = xfp_exp_fx;      
 
        IF( margin_xfp > 0 ) /* xfp_fx was normalized on a larger area than 16ms part of  pcmBuffer  */
        {     
             ASSERT(bfi !=1) ; /* if bfi was set the margin_xfp  does not reflect the correct 16ms part of pcm_buf hist, prev_synth */
             dn_scale =  s_max(0, sub(e_tot_headroom[fs_idx], margin_xfp)); 

             exp_shift = sub(e_tot_headroom[fs_idx], dn_scale);
             used_xfp_exp_fx = sub(xfp_exp_fx, exp_shift); /* the virtual change of the xfp_buffer due to reduced downscaling in L_tot calc  */
        }
        
        /* use semifixed dn_scale as adjusted by margin_xfp in 16 ms region */
        exp_out = xfp_exp_fx;     move16();
        L_tot = winEnCalc(xfp_fx, dn_scale , PhECU_wins[fs_idx][0], rectLengthTab[fs_idx], hamm_len2Tab[fs_idx], &exp_out );
           
        *L_old_xfp_w_E_fx   = L_tot;      move32();

        *old_xfp_w_E_exp_fx = used_xfp_exp_fx   ;  move16();   
       /* this now needs to be in Q1 , used_fx_exp , (exp_out-1-2*e_tot_headroom[fs_idx])/2  */

        *old_Ltot_exp_fx  = exp_out;  /* new proper _Ltot value from winEnCalc function */ 

       
         /* use true Word32 exponent of L_tot */
            

        /* restart oold and old from same  state for init or prevBFI cases  */   
        logic16();
        IF (sub_sat(*oold_xfp_w_E_exp_fx, UNINIT_OR_UNSAFE_OOLD_SENTINEL)  <= 0  ||  /* old xfp_Exp */
            sub_sat(*oold_Ltot_exp_fx, UNINIT_OR_UNSAFE_OOLD_SENTINEL)  <= 0     )    /* new L_tot_exp */
        {
            *L_oold_xfp_w_E_fx   = L_tot;       move32();
            *oold_xfp_w_E_exp_fx = used_xfp_exp_fx;  move16();  
            *oold_Ltot_exp_fx    = *old_Ltot_exp_fx;     /* use   Ltot exp value */ 
        }
    }

    BASOP_sub_sub_end();
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}



