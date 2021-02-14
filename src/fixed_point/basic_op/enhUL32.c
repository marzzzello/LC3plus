/*
  ===========================================================================
   File: ENHUL32.C                                        v.0.6 - 04.Nov.2016
  ===========================================================================

            ENHANCED UNSIGNED 32-BIT ARITHMETIC OPERATORS
   History:
 
  ============================================================================
*/

 /*****************************************************************************
 *
 *  Enhanced Unsigned 32 bit operators :
 *   see complete list in .h file
 * 
 *****************************************************************************/

/*****************************************************************************
 *
 *  Include-Files
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "stl.h"
#include "enhUL32.h"

                       
#if (WMOPS)
extern BASIC_OP multiCounter[MAXCOUNTERS];
extern int currCounter;
#endif /* if WMOPS */


/*****************************************************************************
 *
 *   Constants and Globals
 *
 *****************************************************************************/


/*****************************************************************************
 *
 *   Functions
 *
 *****************************************************************************/
#ifdef ENHUL32

/*___________________________________________________________________________
 |                                                                           |
 |   Function Name : UL_deposit_l                                            |
 |                                                                           |
 |   Purpose :                                                               |
 |                                                                           |
 |   Deposit the 16 bit var1 into the 16 LS bits of the 32 bit output. The   |
 |   16 MS bits of the output are not sign extended.                         |
 |                                                                           |
 |___________________________________________________________________________ */

UWord32 UL_deposit_l(UWord16 uvar){
   UWord32 UL_result;
   UL_result = (UWord32)uvar;   /* no sign extension*/
#if (WMOPS)
    multiCounter[currCounter].L_deposit_l++;
#endif
   return (UL_result);
}


/*___________________________________________________________________________
 |                                                                           |
 |   Function Name : norm_ul                                                  |
 |                                                                           |
 |   Purpose :                                                               |
 |                                                                           |
 |   Produces the number of left shifts needed to normalize the 32 bit varia-|
 |   ble UL_var1 for positive values on the interval with minimum of         |
 |   0 and maximum of 0xffffffff, ; in order to normalize the                | 
 |   result, the following operation must be done :                          |
 |                                                                           |
 |      norm_UL_var1 = UL_lshl(UL_var1, norm_ul(UL_var1)).                    |
 |                                                                           |
 |   Complexity weight : 1                                                   |
 |                                                                           |
 |   Inputs :                                                                |
 |                                                                           |
 |    UL_var1                                                                |
 |          32 bit long unsigned integer (UWord32) whose value falls in the  |
 |          range : 0x0000 0000 <= var1 <= 0xffff ffff.                   |
 |                                                                           |
 |   Outputs :                                                               |
 |                                                                           |
 |    none                                                                   |
 |                                                                           |
 |   Return Value :                                                          |
 |                                                                           |
 |    var_out                                                                |
 |             16 bit short signed integer (Word16) whose value falls in the |
 |             range : 0x0000 0000 <= var_out <= 0x0000 001f. (0..31d)       |
 |___________________________________________________________________________|
*/
Word16 norm_ul( UWord32 UL_var1)
{
    Word16 var_out;

    if (UL_var1 == 0)
    {
        var_out = 0;
    }
    else
    {
            /* simply test shift up until  highest bit is set */
            for (var_out = 0; UL_var1 < (UWord32) 0x80000000U; var_out++)
            {
                UL_var1 <<= 1;
            }       
    }

#if (WMOPS)
    multiCounter[currCounter].norm_l++;  /* now new counter added for UL_ */
#endif
    return (var_out);
}


#ifdef  DRAFT_UL_UPDATE

/*___________________________________________________________________________
 |                                                                           |
 |   Function : UL_addNs                                                     |
 |                                                                           |
 |    addition of  two unsigned 32 bit inputs                                |
 |    no saturation                                                          |
 |                                                                           |
 |    Outputs :                                                              |
 |        *wrap    =  1 if wrap occured, otherwize 0                         |
 |                                                                           |
 |    Return Value :                                                         | 
 |        UL_var3 =  modulo(UL_var1+UL_var2,32)                              | 
 |___________________________________________________________________________|
*/

UWord32  UL_addNs(UWord32 UL_var1, UWord32 UL_var2, UWord16 * wrap )
{
  UWord32 UL_var3; 
  
  /*  no update of Overflow flag  */      
  UL_var3 = UL_var1 + UL_var2;  /* 32 bit wrap  may occur, like in C */
  *wrap = 0; 
  if ( ( ((UWord64)UL_var1 + (UWord64)UL_var2) ) >  ((UWord64) maxUWord32 )  )  
  { 
     *wrap = 1;  /* wrapped output */
  }
#if WMOPS
  multiCounter[currCounter].L_sub++; /* no new counter added */
#endif
  return  UL_var3;
}

#else
/*___________________________________________________________________________
 |                                                                           |
 |   Function Name : UL_addNs                                                |
 |                                                                           |
 |   Purpose :                                                               |
 |                                                                           |
 |   32 bits addition of the two unsigned 32 bits variables                  |
 |   (L_var1+L_var2) with overflow control, but without saturation           |
 |
 |   Outputs :                                                               |
 |        *wrap    =  1 if wrap occured, otherwize 0                                                                      |
 |                                        
 |    Return Value :   
 |        UL_var3 =  modulo(UL_var1+UL_var2,32)  
 |___________________________________________________________________________|
*/

UWord32  UL_addNs(UWord32 UL_var1, UWord32 UL_var2, UWord16* wrap )
{
  UWord32 UL_var3; 
  
  /*  STL Overflow flag is not updated */ 
     
  UL_var3 = UL_var1 + UL_var2;        /* 32bit wrap  may occur, like in C */
   
  if ( ( ((UWord64)UL_var1 + (UWord64)UL_var2) ) >  ((UWord64) maxUWord32 )  )  
  { 
     *wrap = 1;  /* wrapped output */
  }
  else
  {
    *wrap = 0; 
  }

#if WMOPS
  multiCounter[currCounter].L_add++; /* now new counter added for UL_ */
#endif
  return  UL_var3;
}
#endif


#ifdef DRAFT_UL_UPDATE
/*___________________________________________________________________________
 |                                                                           |
 |   Function Name : UL_subNs                                                |
 |                                                                           |                                                          
 |   unsigned  subtraction between two 32 bit values                         |
 |   no saturation                                                           |
 |                                                                           |
 |   Outputs :                                                               |
 |        *wrap    =  1 if wrap (to "negative" occured, otherwise 0          |   
 |                                                                           |
 |    Return Value :                                                         |    
 |        UL_var3 =  modulo(UL_var1-UL_var2, 32)                             | 
 |___________________________________________________________________________|
*/

UWord32 UL_subNs(UWord32 UL_var1, UWord32 UL_var2, UWord16* wrap)
{
  UWord32 UL_var3;

  UL_var3 = UL_var1 - UL_var2;  /*  wrap  may occur, like in C */
  *wrap = 0;                      
  if (UL_var1 < UL_var2)
  { 
    *wrap = 1;  /* "negative" , wrapped output */
  }

#if WMOPS
  multiCounter[currCounter].L_add++;   /* no new counter added  */
#endif
  return  UL_var3;
}
#else
/*___________________________________________________________________________
 |                                                                           |
 |   Function Name : UL_subNs                                                 |
 |                                                                           |
 |   Purpose :                                                               |
 |                                                                           |
 |   32 bits subtraction of the two unsigned 32 bits variables               |
 |   (L_var1-L_var2) with  overflow control, but without saturation          |
 |
 |   Outputs :                                                               |
 |        *sgn    =  1 if wrap (to "negative" occured, otherwise 0                                                                      |
 |                                        
 |    Return Value :   
 |        UL_var3 =  modulo(UL_var1-UL_var2,32) 
 |___________________________________________________________________________|
*/

UWord32 UL_subNs(UWord32 UL_var1, UWord32 UL_var2, UWord16 * sgn)
{
  UWord32 UL_var3;

  UL_var3 = UL_var1 - UL_var2;  /*wrap  may occur, like in C */
  if (UL_var1 >= UL_var2)
  { 
    *sgn = 0;
  }
  else
  {
    *sgn = 1; /* "negative" , wrapped output */
  }

#if WMOPS
  multiCounter[currCounter].L_sub++;  /* now new counter added for UL_ */
#endif
  return  UL_var3;
}

#endif


/*****************************************************************************
 *
 *  Function Name :  Mpy_32_16_uu
 *
 *  Purpose :
 *
 *    Multiplies the 2 unsigned values UL_var1 and uvar2 
 *    with saturation control on 48-bit.( does not happen for unsigned)
 *    The operation is performed in fractional mode :
 *    - UL_var1 is supposed to be in Q32 format.
 *    - var2   is supposed to be in Q16 format.
 *    - The result is produced in Q48 format : UL_varout_h points to the
 *      32 MSBits while varout_l points to the 16 LSBits.
 *
 *  Complexity weight : 2
 *
 *  Inputs :
 *
 *    UL_var1     32 bit long unsigned integer (UWord32) whose value falls in
 *                the range : 0x0000 0000 <= L_var1 <= 0xffff ffff.
 *
 *    var2        16 bit short unsigned integer (UWord16) whose value falls in 
 *                the range : 0x0000 <= var2 <= 0x0000 ffff.
 *
 *  Outputs :
 *
 *    *UL_varout_h 32 bit long unsigned integer (UWord32) whose value falls in
 *                 the range : 0x0000 0000 <= UL_varout_h <= 0xffff ffff.
 *
 *    *varout_l   16 bit short unsigned integer (UWord16) whose value falls in
 *                the range : 0x0000 0000 <= varout_l <= 0x0000 ffff.
 *
 *  Return Value :
 *
 *    none
 *
 *****************************************************************************/
void Mpy_32_16_uu( UWord32 UL_var1, UWord16 uvar2, UWord32 *UL_varout_h, UWord16 *varout_l) {  
   UWord64  UL64_var1; 

   /* 4294967295  * 65535  <  281474976710655 */
   /* (uint64(2)^16-1 )*(uint64(2)^32-1) < (uint64(2)^(16+32)-1) */
   { 
      UL64_var1   = ( UWord64) UL_var1 * ( UWord64)uvar2 ;

      *varout_l   = ( UWord16)( UL64_var1 );

      *UL_varout_h = (UWord32)(UL64_var1>>16);
   }

   #if (WMOPS)
     multiCounter[currCounter].Mpy_32_16_ss++;   /* now new counter added for UL_ */
   #endif /* if WMOPS */

   return;
}


/*****************************************************************************
 *
 *  Function Name :  Mpy_32_32_uu
 *
 *  Purpose :
 *
 *    Multiplies the 2 unsigned values UL_var1 and UL_var2 
 *    with saturation control on 64-bit. (not needed for unsigned) 
 *    The operation is performed in fractional mode :
 *    - UL_var1 and UL_var2 are supposed to be in Q32 format.
 *    - The result is produced in Q64 format : UL_varout_h points to the
 *      32 MSBits while UL_varout_l points to the 32 LSBits.
 *
 *  Complexity weight : 4
 *
 *  Inputs :
 *
 *    UL_var1      32 bit long unsigned integer (UWord32) whose value falls in the
 *                range : 0x0000 0000  <= L_var1 <= 0xffff ffff.
 *
 *    UL_var2      32 bit long unsigned integer (UWord32) whose value falls in the
 *                range : 0x0000 0000  <= L_var2 <= 0xffff ffff.
 *
 *  Outputs :
 *
 *    *UL_varout_h 32 bit long signed integer (Word32) whose value falls in
 *                 the range : 0x0000 0000 <= UL_varout_h <= 0xffff ffff.
 *
 *    *UL_varout_l 32 bit short unsigned integer (UWord32) whose value falls in
 *                 the range : 0x0000 0000 <= UL_varout_l <= 0xffff ffff.
 *
 *
 *  Return Value :
 *
 *    none
 *
 *****************************************************************************/
void Mpy_32_32_uu( UWord32 UL_var1, UWord32 UL_var2, UWord32 *UL_varout_h, UWord32 *UL_varout_l) 
{
    UWord64 UL64_var1;

 /* (uint64(2)^32-1 )*(uint64(2)^32-1) < (uint64(2)^(32+32)-1) */
   {
      UL64_var1     = ((UWord64) UL_var1)*((UWord64) UL_var2); 
      *UL_varout_h  = (UWord32)(UL64_var1>>32);     
      *UL_varout_l  = (UWord32)(UL64_var1);     
   }
   #if (WMOPS)
     multiCounter[currCounter].Mpy_32_32_ss++;   /* now new counter added for UL_ */
   #endif /* if WMOPS */

   return;
}

/*****************************************************************************
 *
 *  Function Name :  UL_Mpy_32_32
 *
 *  Purpose :
 *
 *    Multiplies the 2 unsigned values UL_var1 and UL_var2 
 *    and returns the lower 32 bits, without saturation control. 
 *   
 *    - UL_var1 and UL_var2 are supposed to be in Q32 format.
 *    - The result is produced in Q64 format,  (the 32 LSBits. )
 * 
 *      operates like a  regular 32-by-32 bit unsigned int  multiplication in ANSI-C.
 *          UWord32) = (unsigned int)*(unsigned int);
 *                 
 *
 *  Complexity weight : 2
 *
 *  Inputs :
 *
 *    UL_var1     32 bit long unsigned integer (UWord32) whose value falls in the
 *                range : 0x0000 0000  <= UL_var1 <= 0xffff ffff.
 *
 *    UL_var2     32 bit long unsigned integer (UWord32) whose value falls in the
 *                range : 0x0000 0000  <= UL_var2 <= 0xffff ffff.
 *
 *  Outputs :
 *
 *  Return Value : 
 *     *UL_varout_l 32 bit short unsigned integer (UWord32) whose value falls in
 *      the range : 0x0000 0000 <= UL_varout_l <= 0xffff ffff.
 *
 *    none
 *
 *****************************************************************************/
UWord32  UL_Mpy_32_32( UWord32 UL_var1, UWord32 UL_var2) 
{
    UWord32 UL_varout_l;

#define MASK32  0xFFFFFFFFU
    /* MASK32 may be needed in case Hardware is using larger than 32 bits for UWord32 type) */  
    UL_varout_l       =  (UL_var1&MASK32)*(UL_var2&MASK32);      
    UL_varout_l       =  (UL_varout_l&MASK32);
#undef MASK32      


#if (WMOPS)
    multiCounter[currCounter].Mpy_32_16_ss++;   /* now new counters added for UL_  ops*/
#endif /* if WMOPS */

    return UL_varout_l;
}

#ifdef STL_TYPECASTS 
/*     (Reuse of existing signed STL "L" operators) with
        typecasting  to  make the resulting "UL" code a lot cleaner and more readable. */

UWord32 UL_lshl( UWord32 UL_var1, Word16 var2) {
   return( (UWord32)L_lshl( (Word32) UL_var1, var2));
}

UWord32 UL_lshr( UWord32 UL_var1, Word16 var2) {
   return( (UWord32)L_lshr( (Word32) UL_var1, var2) );
}

UWord32 UL_and(UWord32 UL_var1, UWord32 UL_var2 )
{
   return (UWord32) L_and((Word32)UL_var1,(Word32) UL_var2);
}

UWord32 UL_or(UWord32 UL_var1, UWord32 UL_var2 )
{
   return (UWord32) L_or((Word32)UL_var1,(Word32) UL_var2);
}

UWord32 UL_xor(UWord32 UL_var1, UWord32 UL_var2 )
{
   return (UWord32) L_xor((Word32)UL_var1,(Word32) UL_var2);
}

UWord32 UL_deposit_h(UWord16 uvar1)
{
   return (UWord32) L_deposit_h((Word32)uvar1);
}

UWord16 u_extract_h(UWord32 UL_var1)
{
   return (UWord16) extract_h((Word32)UL_var1);
}

UWord16 u_extract_l(UWord32 UL_var1)
{
   return (UWord32)extract_l((Word32)UL_var1);
}

/* enable convenient reuse of Non-saturating UL_subNs , UL_addNs  while "D"iscarding the sgn/wrap output flags */ 
UWord32 UL_subNsD(UWord32 UL_var1, UWord32 UL_var2 )
{ 
    UWord16 dummy_sign;
    return UL_subNs(UL_var1,UL_var2,&dummy_sign );
}

UWord32 UL_addNsD(UWord32 UL_var1, UWord32 UL_var2 )
{  
   UWord16 dummy_wrap;
   return UL_addNs(UL_var1,UL_var2,&dummy_wrap );
}
#endif

#endif /* ENHUL32 */

/* end of file */
