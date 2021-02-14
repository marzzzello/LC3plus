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
#include "stl.h"
#include "basop32.h"



void processPLCupdate_fx(AplcSetup *plcAd, Word16 x_fx[], Word16 q_fx_exp, Word16 concealMethod, Word16 frame_length,
                         Word16 fs_idx, Word16 *nbLostFramesInRow, Word16 *prev_prev_bfi, Word16 *prev_bfi, Word16 bfi, Word16 scf_q[],
                         Word16 ola_mem_fx[], Word16 ola_mem_fx_exp, Word16 *ns_cum_alpha)
{

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("process_plc_update_fx",         
        0
                   );
#endif


    processPLCUpdateAfterIMDCT_fx(x_fx, q_fx_exp, concealMethod, frame_length, fs_idx, nbLostFramesInRow, prev_prev_bfi, prev_bfi, bfi,
                                  scf_q, ns_cum_alpha, plcAd); /* NB *prev_bfi updated here */

    UNUSED(ola_mem_fx);
    UNUSED(ola_mem_fx_exp);
    IF ( plcAd != 0 )
    { 
        /*  reuse/inplace the most recent 16 ms of x_old_tot without additional rescaling,  keep exponent aligned with tdc pitch buffer to  save WMOPS */
        ASSERT( (&plcAd->x_old_tot_fx[plcAd->max_len_pcm_plc - LprotSzPtr[fs_idx] ])  == plcAd->PhECU_xfp_fx );   
        plcAd->PhECU_xfp_exp_fx  = plcAd->q_fx_old_exp;    move16(); /* exponent used by concealmethod 2 in prevBfi frames and also right after  non bfi frames */
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif

}

void processPLCupdateSpec_fx(Word16 q_old_d_fx[], Word16 *q_old_fx_exp, Word32 q_d_fx[], Word16 *q_fx_exp, Word16 yLen)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
        Word16  s;
    );

    /* save spectrum and the corresponding exponent */
    s             = getScaleFactor32(q_d_fx, yLen);

    *q_old_fx_exp = sub(*q_fx_exp, s);


    FOR (i = 0; i < yLen; i++)
    {
        q_old_d_fx[i] = round_fx_sat(L_shl_sat(q_d_fx[i], s)); /*  */
    }

    Dyn_Mem_Deluxe_Out();
}



void processPLCspec2shape_fx(Word16 prev_bfi, Word16 bfi, Word16 q_old_d_fx[], Word16 yLen,  
                             Word16 *stPhECU_oold_grp_shape_fx, Word16 *stPhECU_old_grp_shape_fx)

#define L_GRP_DC 4 

{
    Counter i,l; 
    Word16  *pX, tmp; 
    Word16  N_grp,l_grp;
    Word16  man, expo;
    Word32  L_acc;
    Word32  L_tot;
    Word32  L_grp_shape[MAX_LGW];
    Word16  grp_shape[MAX_LGW]; /**/
    Word16  fs_idx,local_prev_bfi;
   

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("process_plc_spec_2_shape_fx", sizeof(struct {
                   Counter i,l;
                   Word16 *pX; 
                   Word16  N_grp,l_grp;
                   Word32  L_acc;
                   Word32  L_tot;
                   Word32  L_grp_shape[MAX_LGW];
                   Word16  fs_idx,local_prev_bfi;
               }));
#endif
    BASOP_sub_sub_start("PhECU::GF::process_plc_spec_2_shape_fx");

    IF(sub(bfi, 1) != 0)
    {
       fs_idx = mult(yLen, (Word16)(32768.0 / (99.0))); /* truncation needed , i.e no rounding can be applied here */
       N_grp = xavg_N_grp_fx[fs_idx];   move16();

       local_prev_bfi = prev_bfi; move16();
       if (sub(local_prev_bfi, 2)==0) 
       {
          local_prev_bfi = 0; move16();
       }

       if( stPhECU_oold_grp_shape_fx[0] < 0 )  
       {
          local_prev_bfi = 1 ;   move16(); /* handle startup in the case 2nd frame is a  BFI frame */     
       }

        /* Copy old to oold grp shape */
        FOR( i=0; i < MAX_LGW ; i++) 
        {
            stPhECU_oold_grp_shape_fx[i] = stPhECU_old_grp_shape_fx[i];  move16();
        }     
   
 
        /* Accumulate DC bin(s) to total */
        pX    = q_old_d_fx;          /*  ptr setup */ 
        L_tot = L_deposit_h(0);      /* acc on negative side */
       
        FOR( i= 0; i < L_GRP_DC; i++)
        {
            tmp    = shr_pos(*pX++ ,spec_shape_headroom[fs_idx]); /* scale down upscaled MDCT to create some headroom */
            L_tot  = L_msu0(L_tot, tmp, tmp);       
            
        }
 
        /* Accumulate middle subbands and add to total */
        FOR( i=0; i < sub(N_grp,1) ; i++)
        {
             
            L_acc  = L_deposit_h(0);  /* acc on negative side */
            l_grp  = sub(mdct_grp_bins_fx[i+1], mdct_grp_bins_fx[i]);   move16();

         

            FOR(l=0;l<l_grp; l++)
            {                           
                tmp    = shr(*pX++ ,spec_shape_headroom[fs_idx]);   
                L_acc  = L_msu0(L_acc, tmp, tmp);  
               
            }
      
 
            L_grp_shape[i] = L_negate(L_acc);                   move32();
            L_tot          = L_add(L_tot, L_acc);               move32();   /* two negative numbers added */   
        }

        /* Accumulate last subbband and add to total */
 
        L_acc = L_deposit_h(0);    
        l_grp = sub(sub(mdct_grp_bins_fx[N_grp], mdct_grp_bins_fx[N_grp-1]),L_GRP_DC);  

 
        FOR(l=0; l<l_grp; l++)
        {
            tmp   = shr(*pX++, spec_shape_headroom[fs_idx]); 
            L_acc = L_msu0(L_acc, tmp, tmp);
            
        }

 
        L_grp_shape[sub(N_grp,1)] = L_negate(L_acc);   move32();


        L_tot                = L_add(L_tot, L_acc);       /* two negative numbers added */   
       
        L_tot  = L_max( -(INT_MAX), L_tot); /* conditionally add 1 to negative side, to avoid possible saturation in L_negate */ 
        L_tot  = L_negate(L_tot);           /* no saturation here as L_tot is != INT_MIN */
                   

        /* Normalize shape */
        /* norm_scale = 1/L_tot; */

        IF (L_tot > 0)
        {
            FOR(i=0; i < N_grp ; i++)
            {
                man          = plc_phEcu_ratio_fx(L_grp_shape[i], L_tot, &expo); /* The mantissa is considered in Q15 output in Q14 */
                grp_shape[i] = shr_sat(man, sub(expo,1));    /* gfrom Q14 to in Q15 (Due to saturation, it is automatically bound inside [-1.0,1.0].) */
            }
        }
        ELSE
        {
            FOR(i=0; i < N_grp ; i++)
            {
                grp_shape[i] = GRP_SHAPE_INIT;                      move16();    
              
            }             
        }

        /* copy to output */
        FOR(i=0; i < N_grp ; i++)  
        {
            stPhECU_old_grp_shape_fx[i] = grp_shape[i];    move16();
        }
        FOR(i = N_grp; i < MAX_LGW ; i++)  
        {
           stPhECU_old_grp_shape_fx[i] = GRP_SHAPE_INIT;    move16();        
        }
 


        /* handle oold states for the frame sequence    BAD, GOOD,  NEXT_BAD */
        if(sub(local_prev_bfi, 1)==0)
        {
            FOR( i=0; i < MAX_LGW ; i++)   
            {
                stPhECU_oold_grp_shape_fx[i] = stPhECU_old_grp_shape_fx[i] ;  move16();
            }
        }

    }

    BASOP_sub_sub_end();
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif

}




