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

/* initilize a short vector */
void plc_phEcu_initWord16(Word16 *     vec,   /*i/o : short vector pointer       */
                          const Word16 value, /*i   : short initialization value */
                          const Word16 len)   /*i   : number of elements         */
{
    Counter n;

    FOR (n = 0; n < len; n++)
    {
        vec[n] = value; move16();
    }
}

/* scale inplace with allowed saturation in upscaling ,  function not available in basop_util */
void Scale_sig_sat(Word16       x[], /* i/o: signal to scale                 Qx        */
                  const Word16 lg,   /* i  : size of x[]                     Q0        */
                  const Word16 exp0  /* i  : exponent: x = round(x << exp)   Qx ?exp  */
)
{
    Counter i;
    Word16  tmp;
    IF (exp0 > 0)
    {
        FOR (i = 0; i < lg; i++)
        {
            x[i] = shl_sat(x[i], exp0);  move16(); /* no saturation warnings triggered here */
        }
        return;
    }
    IF (exp0 < 0)
    {
        tmp = shl(-32768, s_max(exp0, -15)); /* we use negative to correctly represent 1.0 */
        FOR (i = 0; i < lg; i++)
        {
            x[i] = msu_r(0, x[i], tmp);  move16(); /* msu instead of mac because factor is negative */
        }
        return;
    }
}

void plc_phEcu_minval_fx(const Word16 *inp,      /* i  : vector       */
                         const Word16  len,      /* i  : length       */
                         Word16       *minvalPtr /* o  : min  value Ptr    */
)
{
    Word16  minTmp;
    Counter pos;

    minTmp = inp[0];  move16();
    assert(len>1);
    FOR (pos = 1; pos < len; pos++)
    {
        minTmp = s_min(inp[pos], minTmp);
    }

    *minvalPtr = minTmp;  move16();
}

void plc_phEcu_maxval_fx(const Word16 *inp,      /* i  : vector     */
                         const Word16  len,      /* i  : length     */
                         Word16       *maxvalPtr /* o  : *maxvalPtr */
)
{
    Word16  maxTmp;
    Counter pos;

    maxTmp = inp[0];  move16();

    assert(len>1);
    FOR (pos = 1; pos < len; pos++)
    {
        maxTmp = s_max(inp[pos], maxTmp);
    }
    *maxvalPtr = maxTmp;  move16();
}

/* in case a  value (e.g max or min)  is already known , find the first corresponding array  index */
Word16  plc_phEcu_find_ind_fx(                        /* o  : output maximum  indx 0.. len-1    */
                              const Word16 *inp,      /* i  : vector     */
                              const Word16  len,      /* i  : length     */
                              const Word16  val   /* i  : value to find     */
)       
{
   Word16  val_ind;
   Counter pos;

   val_ind = -1;  move16();

   FOR(pos = 0; pos < len; pos++)
   {
      if (sub(inp[pos], val) == 0)
      {
         val_ind = pos;  move16();
      }
   }
   return   val_ind;
}

/*-----------------------------------------------------------------------------
 * ratio_fx()
 *
 * Divide the numerator by the denominator.
 *----------------------------------------------------------------------------*/
Word16 plc_phEcu_ratio_fx(                     /* o : quotient   in Q14       */
                          const Word32 numer,  /* i : numerator               */
                          const Word32 denom,  /* i : denominator             */
                          Word16 *expo)        /* o : req shift of quotient   */
{
    Word16 expNumer, expDenom;
    Word16 manNumer, manDenom;
    Word16 quotient;
#ifdef DYNMEM_COUNT
    Dyn_Mem_In("plc_phEcu_ratio_fx", sizeof(struct {
                   Word16 expNumer, expDenom;
                   Word16 manNumer, manDenom;
                   Word16 quotient;
               }));
#endif

    expDenom = norm_l(denom);                     /* exponent */
    manDenom = extract_h(L_shl(denom, expDenom)); /* mantissa */
    expNumer = norm_l(numer);                     /* exponent */
    manNumer = extract_h(L_shl(numer, expNumer)); /* mantissa */
    manNumer = shr_pos(manNumer, 1);              /* Ensure the numerator < the denominator */
    quotient = div_s(manNumer, manDenom);         /* in Q14 */

    *expo = sub(expNumer, expDenom);
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
    return quotient; /* Q14 */
}

Word32 winEnCalc(                          /* o:  output summed energy Ltot */
                 const Word16 *x,          /* i: Input signal */
                 const Word16 headroom_shift,    /* i: headroom_shift */
                 const Word16 *win,        /* i: left side Window coefficients */
                 const Word16  rectLength, /* i: Offset in between the 1st and 2nd symmetric halves of the Hamming window */
                 const Word16  halfLength,  /* i: Half of the total length of a complete Hamming window. */
                 Word16    *exp             /* i/o : i exp of Word16 variable x ,  o:Lexp of output Word32 sum */     
                 )
{
    Counter       i;
    Word32        L_tot; 
    const Word16 *pX, *pW;
    Word16        tmp, tmp_RL;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("PhEcu::GF::winEnCalc", sizeof(struct {
                   Counter       i;
                    Word32        L_tot;               
                   const Word16 *pX, *pW;
                   Word16 tmp, tmp_RL;
               }));
#endif


    BASOP_sub_sub_start("PhECU::GF::winEnCalc");

    L_tot = INT_MAX; move32(); /*acc is on negative side , but as all accumulatio is positive, we make use of one extra bit   */
    pX   = x;
    pW   = win;
   

    assert( headroom_shift>=0 );
    FOR (i = 0; i < halfLength; i++) /* 1st symmetric half of the Hamming window */
    {
          tmp   = mult(*pX++, *pW++);   
          tmp   = shr_pos(tmp, headroom_shift);  /* shr  may/create  bias on the negative side , costly options are shr_r or use msu_r */
          L_tot = L_msu0(L_tot, tmp, tmp);       /* acc on negative energy side */
    }

    /* Periodic filter - one more rect sample before end tapering */
    tmp_RL = add(rectLength, 1);
    ASSERT(rectLength != 0);

    FOR (i = 0; i < tmp_RL; i++) /* If rectLength is zero, it's a pure Hamming window; otherwise Hamming-Rectangular. */
    {    
          tmp   = shr_pos( *pX++, headroom_shift);
          L_tot = L_msu0(L_tot, tmp, tmp); /* acc on negative side */
    }

    tmp_RL = sub(halfLength, 1);
    ASSERT(rectLength != 0);

    FOR (i = 0; i < tmp_RL; i++) /* 2nd symmetric half of the Hamming window. */
    {
        tmp   = mult(*pX++, *(--pW));   
        tmp   = shr_pos(tmp, headroom_shift);
        L_tot = L_msu0(L_tot, tmp, tmp);
    }

  /*  Lexp = 2*(incoming_exp + dnshift) + 1  , 2x for square + 1(for msu0 DSP dn shift)*/
   *exp   =   add(shl_pos(add(*exp, headroom_shift),1),1);  

     /* handle wrap on zero point */
    IF( L_tot >= 0 )
    {  /* L_tot positive           --> less than 32 bits needed, */
       L_tot = L_add(L_tot,(INT_MIN+1));
       if( L_tot == 0 )
       {      
          *exp =  LTOT_MIN_EXP; /* WC is actually (-(15+4)*2 + 1 +1  -31) */ ;  move16();
       }  
       L_tot =  L_min(L_tot, -1);     /* make sure there is energy for future ratio calculations */
    } 
    ELSE 
    {   /* L_tot negative --> more than 31 bits needed for sum , scale  32 bit sum within 31 bits  and  adjust exp */
 
         L_tot = L_shr_pos(L_add(L_tot,1),1);  /* rnd by adding 1,  then use 50% contribution from negative  side */
         L_tot = L_add(L_tot, INT_MIN>>1);     /* add 50% contribution from positive side */      
    
        *exp =  add(*exp, 1);             move16(); 
    } 

    L_tot = L_max( -(INT_MAX), L_tot);  /* guard against max accumulation on the  negative side , should only occur for rectangle windows */
    L_tot = L_negate(L_tot); /* no saturation here */

    /* activate when xfp_exp is not used any longer */  
    /*  pre-maximize the mantissa for  the following  steps  in burst_ana_dx  */
    tmp   = norm_l(L_tot);
    L_tot =  L_shl(L_tot,tmp);
    *exp  =  sub(*exp, tmp);  move16();

 BASOP_sub_sub_end();
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
   
    return L_tot;
}

