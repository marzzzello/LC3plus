/*
  ===========================================================================
   File: ENH40.C                                         v.2.3 - 30.Nov.2009
  ===========================================================================

            ITU-T  STL  BASIC OPERATORS

            40-BIT ARITHMETIC OPERATORS

   History:
   07 Nov 04   v2.0     Incorporation of new 32-bit / 40-bit / control
                        operators for the ITU-T Standard Tool Library as 
                        described in Geneva, 20-30 January 2004 WP 3/16 Q10/16
                        TD 11 document and subsequent discussions on the
                        wp3audio@yahoogroups.com email reflector.

   31 Mar 15   v2.1E    Removal of operators not used in the EVS codec.

  ============================================================================
*/


/*****************************************************************************
 *
 *  Enhanced 40 bit operators :
 *
 *    Mpy_32_16_ss()
 *    Mpy_32_32_ss()
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

#if (WMOPS)
extern BASIC_OP multiCounter[MAXCOUNTERS];
extern int currCounter;
#endif /* if WMOPS */


/*****************************************************************************
 *
 *  Local Functions
 *
 *****************************************************************************/
static __inline Word40 L40_shr(   Word40 L40_var1, Word16 var2);
static __inline Word40 L40_shl(   Word40 L40_var1, Word16 var2);
static __inline Word40 L40_set( Word40 L40_var1);
static __inline UWord16 Extract40_L( Word40 L40_var1);
static __inline Word40 L40_mult( Word16 var1, Word16 var2);
static __inline Word40 L40_add( Word40 L40_var1, Word40 L40_var2);
static __inline Word40 L40_mac( Word40 L40_var1, Word16 var2, Word16 var3);
static __inline UWord32 L_Extract40( Word40 L40_var1) ;

/*****************************************************************************
 *
 *  Macros for 40 bit arithmetic overflow management :
 *  Upon 40-bit overflow beyond MAX_40 or underflow beyond MIN_40,
 *  the application will exit.
 *
 *****************************************************************************/
#ifndef HIDE_UNUSED_BASOP
#define L40_OVERFLOW_OCCURED(  L40_var1) (Overflow = 1, exit(1), L40_var1)
#define L40_UNDERFLOW_OCCURED( L40_var1) (Overflow = 1, exit(2), L40_var1)
#else
#define L40_OVERFLOW_OCCURED(  L40_var1) (exit(1), L40_var1)
#define L40_UNDERFLOW_OCCURED( L40_var1) (exit(2), L40_var1)
#endif

/*****************************************************************************
 *
 *  Constants and Globals
 *
 *****************************************************************************/


/*****************************************************************************
 *
 *  Functions
 *
 *****************************************************************************/


/*****************************************************************************
 *
 *  Function Name :  Mpy_32_16_ss
 *
 *  Purpose :
 *
 *    Multiplies the 2 signed values L_var1 and var2 with saturation control
 *    on 48-bit. The operation is performed in fractional mode :
 *    - L_var1 is supposed to be in 1Q31 format.
 *    - var2   is supposed to be in 1Q15 format.
 *    - The result is produced in 1Q47 format : L_varout_h points to the
 *      32 MSBits while varout_l points to the 16 LSBits.
 *
 *  Complexity weight : 2
 *
 *  Inputs :
 *
 *    L_var1      32 bit long signed integer (Word32) whose value falls in
 *                the range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
 *
 *    var2        16 bit short signed integer (Word16) whose value falls in
 *                the range : 0xffff 8000 <= var2 <= 0x0000 7fff.
 *
 *  Outputs :
 *
 *    *L_varout_h 32 bit long signed integer (Word32) whose value falls in
 *                the range : 0x8000 0000 <= L_varout_h <= 0x7fff ffff.
 *
 *    *varout_l   16 bit short unsigned integer (UWord16) whose value falls in
 *                the range : 0x0000 0000 <= varout_l <= 0x0000 ffff.
 *
 *  Return Value :
 *
 *    none
 *
 *****************************************************************************/
void Mpy_32_16_ss( Word32 L_var1, Word16 var2, Word32 *L_varout_h, UWord16 *varout_l) {
   Word16 var1_h;
   UWord16 uvar1_l;
   Word40 L40_var1;

   if( (L_var1 == ( Word32) 0x80000000)
    && (var2   == ( Word16) 0x8000)) {
      *L_varout_h = 0x7fffffff;
      *varout_l = ( UWord16) 0xffff;

    } else {
      uvar1_l = extract_l( L_var1);
      var1_h = extract_h( L_var1);

      /* Below line can not overflow, so we can use << instead of L40_shl.  */
      L40_var1 = (( Word40) (( Word32) var2 * ( Word32) uvar1_l)) << 1;

      *varout_l = Extract40_L( L40_var1);

      L40_var1 = L40_shr( L40_var1, 16);
      L40_var1 = L40_mac( L40_var1, var2, var1_h);

      *L_varout_h = L_Extract40( L40_var1);

      #if (WMOPS)
      multiCounter[currCounter].extract_l--;
      multiCounter[currCounter].extract_h--;
      multiCounter[currCounter].Extract40_L--;
      multiCounter[currCounter].L40_shr--;
      multiCounter[currCounter].L40_mac--;
      multiCounter[currCounter].L_Extract40--;
      #endif /* if WMOPS */
   }

   #if (WMOPS)
   multiCounter[currCounter].Mpy_32_16_ss++;
   #endif /* if WMOPS */

   return;
}


/*****************************************************************************
 *
 *  Function Name :  Mpy_32_32_ss
 *
 *  Purpose :
 *
 *    Multiplies the 2 signed values L_var1 and L_var2 with saturation control
 *    on 64-bit. The operation is performed in fractional mode :
 *    - L_var1 and L_var2 are supposed to be in 1Q31 format.
 *    - The result is produced in 1Q63 format : L_varout_h points to the
 *      32 MSBits while L_varout_l points to the 32 LSBits.
 *
 *  Complexity weight : 4
 *
 *  Inputs :
 *
 *    L_var1      32 bit long signed integer (Word32) whose value falls in the
 *                range : 0x8000 0000  <= L_var1 <= 0x7fff ffff.
 *
 *    L_var2      32 bit long signed integer (Word32) whose value falls in the
 *                range : 0x8000 0000  <= L_var2 <= 0x7fff ffff.
 *
 *  Outputs :
 *
 *    *L_varout_h 32 bit long signed integer (Word32) whose value falls in
 *                the range : 0x8000 0000 <= L_varout_h <= 0x7fff ffff.
 *
 *    *L_varout_l 32 bit short unsigned integer (UWord32) whose value falls in
 *                the range : 0x0000 0000 <= L_varout_l <= 0xffff ffff.
 *
 *
 *  Return Value :
 *
 *    none
 *
 *****************************************************************************/
void Mpy_32_32_ss( Word32 L_var1, Word32 L_var2, Word32 *L_varout_h, UWord32 *L_varout_l) {
   UWord16 uvar1_l, uvar2_l;
   Word16   var1_h,  var2_h;
   Word40 L40_var1;

   if( (L_var1 == ( Word32)0x80000000)
    && (L_var2 == ( Word32)0x80000000)) {
      *L_varout_h = 0x7fffffff;
      *L_varout_l = ( UWord32)0xffffffff;

   } else {

      uvar1_l = extract_l( L_var1);
      var1_h = extract_h( L_var1);
      uvar2_l = extract_l( L_var2);
      var2_h = extract_h( L_var2);

      /* Below line can not overflow, so we can use << instead of L40_shl.  */
      L40_var1 = (( Word40) (( UWord32) uvar2_l * ( UWord32) uvar1_l)) << 1;

      *L_varout_l = 0x0000ffff & L_Extract40( L40_var1);

      L40_var1 = L40_shr( L40_var1, 16);
      L40_var1 = L40_add( L40_var1, (( Word40) (( Word32) var2_h * ( Word32) uvar1_l)) << 1);
      L40_var1 = L40_add( L40_var1, (( Word40) (( Word32) var1_h * ( Word32) uvar2_l)) << 1);
      *L_varout_l |= (L_Extract40( L40_var1)) << 16;

      L40_var1 = L40_shr( L40_var1, 16);
      L40_var1 = L40_mac( L40_var1, var1_h, var2_h);

      *L_varout_h = L_Extract40( L40_var1);

      #if (WMOPS)
      multiCounter[currCounter].extract_l-=2;
      multiCounter[currCounter].extract_h-=2;
      multiCounter[currCounter].L_Extract40-=3;
      multiCounter[currCounter].L40_shr-=2;
      multiCounter[currCounter].L40_add-=2;
      multiCounter[currCounter].L40_mac--;
      #endif /* if WMOPS */
   }

   #if (WMOPS)
    multiCounter[currCounter].Mpy_32_32_ss++;
   #endif /* if WMOPS */

   return;
}


/*****************************************************************************
 *
 *  Function Name : L40_set
 *
 *  Purpose :
 *
 *    Assigns a 40 constant to a Word40 with adequate initialization depending
 *    on underlying architecture constraints (for example to keep consistency
 *    of sign bits). Current implementation only validated on MSVC++6.0.
 *
 *  Complexity weight : 3
 *
 *  Inputs :
 *
 *    L40_var1    40 bit long signed integer (Word40) whose value falls in the
 *                range : MIN_40 <= L40_var1 <= MAX_40.
 *
 *  Outputs :
 *
 *    none
 *
 *  Return Value :
 *
 *    L40_var_out 40 bit long signed integer (Word40) whose value falls in
 *                the range : MIN_40 <= L40_var_out <= MAX_40.
 *
 *****************************************************************************/
/*#ifdef _MSC_VER*/
static __inline Word40 L40_set( Word40 L40_var1) {
   Word40 L40_var_out;

#if defined(_MSC_VER) && (_MSC_VER <= 1200)
   L40_var_out =  L40_var1 & 0x000000ffffffffff;

   if( L40_var1 & 0x8000000000)
      L40_var_out = L40_var_out | 0xffffff0000000000;
#else
   L40_var_out =  L40_var1 & 0x000000ffffffffffLL;

   if( L40_var1 & 0x8000000000LL)
      L40_var_out = L40_var_out | 0xffffff0000000000LL;
#endif

   #if WMOPS
   multiCounter[currCounter].L40_set++;
   #endif /* if WMOPS */

   return( L40_var_out);
}
/*#endif*/ /* ifdef _MSC_VER */


/*****************************************************************************
 *
 *  Function Name : Extract40_L
 *
 *  Purpose :
 *
 *    Returns the bits [15-0] of L40_var1.
 *
 *  Complexity weight : 1
 *
 *  Inputs :
 *
 *    L40_var1    40 bit long signed integer (Word40) whose value falls in the
 *                range : MIN_40 <= L40_var1 <= MAX_40.
 *
 *  Outputs :
 *
 *    none
 *
 *  Return Value :
 *
 *    var_out     16 bit short unsigned integer (UWord16) whose value falls in
 *                the range : MIN_U_16 <= var_out <= MAX_U_16.
 *
 *****************************************************************************/
static __inline UWord16 Extract40_L( Word40 L40_var1) {
   UWord16 var_out;

   var_out = ( UWord16)( L40_var1);

   #if (WMOPS)
   multiCounter[currCounter].Extract40_L++;
   #endif /* if WMOPS */

   return( var_out);
}


/*****************************************************************************
 *
 *  Function Name : L40_mult
 *
 *  Purpose :
 *
 *    Multiplies var1 by var2 and shifts the result left by 1. Returns the
 *    full precision result on 40-bit.
 *    L40_mult( var1, var2) = shiftleft(( var1 times var2), 1)
 *
 *  Complexity weight : 1
 *
 *  Inputs :
 *
 *    var1        16 bit short signed integer (Word16) whose value falls in
 *                the range : MIN_16 <= var1 <= MAX_16.
 *
 *    var2        16 bit short signed integer (Word16) whose value falls in
 *                the range : MIN_16 <= var2 <= MAX_16.
 *
 *  Outputs :
 *
 *    none
 *
 *  Return Value :
 *
 *    L40_var_out 40 bit long signed integer (Word40) whose value falls in
 *                the range : MIN_40 <= L40_var_out <= MAX_40.
 *
 *****************************************************************************/
static __inline Word40 L40_mult( Word16 var1, Word16 var2) {
   Word32 L_var_out;
   Word40 L40_var_out;

   L_var_out = ( Word32) var1 * ( Word32) var2;
   L40_var_out = ( Word40) L_var_out;

   /* Below line can not overflow, so we can use << instead of L40_shl. */
   L40_var_out = L40_var_out << 1;

   #if (WMOPS)
   multiCounter[currCounter].L40_mult++;
   #endif /* if WMOPS */

   return( L40_var_out);
}

/*****************************************************************************
 *
 *  Function Name : L40_add
 *
 *  Purpose :
 *
 *    Adds L40_var1 and L40_var2 and returns the 40-bit result.
 *    Calls the macro L40_UNDERFLOW_OCCURED() in case of underflow on 40-bit.
 *    Calls the macro L40_OVERFLOW_OCCURED()  in case of overflow  on 40-bit.
 *
 *  Complexity weight : 1
 *
 *  Inputs :
 *
 *    L40_var1    40 bit long signed integer (Word40) whose value falls in the
 *                range : MIN_40 <= L40_var1 <= MAX_40.
 *
 *    L40_var2    40 bit long signed integer (Word40) whose value falls in the
 *                range : MIN_40 <= L40_var2 <= MAX_40.
 *
 *  Outputs :
 *
 *    none
 *
 *  Return Value :
 *
 *    L40_var_out 40 bit long signed integer (Word40) whose value falls in
 *                the range : MIN_40 <= L40_var_out <= MAX_40.
 *
 *****************************************************************************/
static __inline Word40 L40_add( Word40 L40_var1, Word40 L40_var2) {
   Word40 L40_var_out;

   L40_var_out = L40_var1 + L40_var2;

#if defined(_MSC_VER) && (_MSC_VER <= 1200)
   if(  ((( L40_var1    & 0x8000000000) >> 39) != 0)
     && ((( L40_var2    & 0x8000000000) >> 39) != 0)
     && ((( L40_var_out & 0x8000000000) >> 39) == 0)) {
      L40_var_out = L40_UNDERFLOW_OCCURED( L40_var_out);

   } else if( (((L40_var1    & 0x8000000000) >> 39) == 0)
           && (((L40_var2    & 0x8000000000) >> 39) == 0)
           && (((L40_var_out & 0x8000000000) >> 39) != 0)) {
      L40_var_out = L40_OVERFLOW_OCCURED( L40_var_out);
   }
#else
   if(  ((( L40_var1    & 0x8000000000LL) >> 39) != 0)
     && ((( L40_var2    & 0x8000000000LL) >> 39) != 0)
     && ((( L40_var_out & 0x8000000000LL) >> 39) == 0)) {
      L40_var_out = L40_UNDERFLOW_OCCURED( L40_var_out);

   } else if( (((L40_var1    & 0x8000000000LL) >> 39) == 0)
           && (((L40_var2    & 0x8000000000LL) >> 39) == 0)
           && (((L40_var_out & 0x8000000000LL) >> 39) != 0)) {
      L40_var_out = L40_OVERFLOW_OCCURED( L40_var_out);
   }
#endif

   #if (WMOPS)
   multiCounter[currCounter].L40_add++;
   #endif /* if WMOPS */

   BASOP_CHECK(0, 0);

   return( L40_var_out);
}

/*****************************************************************************
 *
 *  Function Name : L40_mac
 *
 *  Purpose :
 *
 *    Multiplies var2 by var3. Shifts left the 40-bit result by 1 and adds
 *    the result to L40_var1. Returns a 40 bit result.
 *    L40_mac( L40_var1, var2, var3)
 *    = L40_add( L40_var1, L40_mult( var2, var3))
 *    Calls the macro L40_UNDERFLOW_OCCURED() in case of underflow on 40-bit.
 *    Calls the macro L40_OVERFLOW_OCCURED()  in case of overflow  on 40-bit.
 *
 *  Complexity weight : 1
 *
 *  Inputs :
 *
 *    L40_var1    40 bit long signed integer (Word40) whose value falls in the
 *                range : MIN_40 <= L40_var1 <= MAX_40.
 *
 *    var2        16 bit short signed integer (Word16) whose value falls in
 *                the range : MIN_16 <= var2 <= MAX_16.
 *
 *    var3        16 bit short signed integer (Word16) whose value falls in
 *                the range : MIN_16 <= var3 <= MAX_16.
 *
 *  Outputs :
 *
 *    none
 *
 *  Return Value :
 *
 *    L40_var_out 40 bit long signed integer (Word40) whose value falls in
 *                the range : MIN_40 <= L40_var_out <= MAX_40.
 *
 *****************************************************************************/
static __inline Word40 L40_mac( Word40 L40_var1, Word16 var2, Word16 var3) {
   Word40 L40_var_out;

   L40_var_out = L40_mult( var2, var3);
   L40_var_out = L40_add( L40_var1, L40_var_out);

   #if (WMOPS)
   multiCounter[currCounter].L40_mult--;
   multiCounter[currCounter].L40_add--;
   multiCounter[currCounter].L40_mac++;
   #endif /* if WMOPS */

   return( L40_var_out);
}

/*****************************************************************************
 *
 *  Function Name : L_Extract40
 *
 *  Purpose :
 *
 *    Returns the bits [31-0] of L40_var1.
 *
 *  Complexity weight : 1
 *
 *  Inputs :
 *
 *    L40_var1    40 bit long signed integer (Word40) whose value falls in the
 *                range : MIN_40 <= L40_var1 <= MAX_40.
 *
 *  Outputs :
 *
 *    none
 *
 *  Return Value :
 *
 *    L_var_out   32 bit long unsigned integer (UWord32) whose value falls in
 *                range : MIN_U_32 <= L_var_out <= MAX_U_32.
 *
 *****************************************************************************/
static __inline UWord32 L_Extract40( Word40 L40_var1) {
   UWord32 L_var_out;

   L_var_out = ( UWord32) L40_var1;

   #if (WMOPS)
   multiCounter[currCounter].L_Extract40++;
   #endif /* if WMOPS */

   return(L_var_out);
}


/*****************************************************************************
 *
 *  Function Name : L40_shl
 *
 *  Purpose :
 *
 *    Arithmetically shifts left L40_var1 by var2 positions.
 *    - If var2 is negative, L40_var1 is shifted to the LSBits by (-var2)
 *      positions with extension of the sign bit.
 *    - If var2 is positive, L40_var1 is shifted to the MSBits by (var2)
 *      positions.
 *    Calls the macro L40_UNDERFLOW_OCCURED() in case of underflow on 40-bit.
 *    Calls the macro L40_OVERFLOW_OCCURED()  in case of overflow  on 40-bit.
 *
 *  Complexity weight : 1
 *
 *  Inputs :
 *
 *    L40_var1    40 bit long signed integer (Word40) whose value falls in the
 *                range : MIN_40 <= L40_var1 <= MAX_40.
 *
 *    var2        16 bit short signed integer (Word16) whose value falls in
 *                the range : MIN_16 <= var2 <= MAX_16.
 *
 *  Outputs :
 *
 *    none
 *
 *  Return Value :
 *
 *    L40_var_out 40 bit long signed integer (Word40) whose value falls in
 *                the range : MIN_40 <= L40_var_out <= MAX_40.
 *
 *****************************************************************************/
static __inline Word40 L40_shl( Word40 L40_var1, Word16 var2) {

   Word40 L40_var_out;
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
   Word40 L40_constant = L40_set( 0xc000000000);
#else
   Word40 L40_constant = L40_set( 0xc000000000LL);
#endif

   if( var2 < 0) {
      var2 = -var2;
      L40_var_out = L40_shr( L40_var1, var2);

      #if (WMOPS)
      multiCounter[currCounter].L40_shr--;
      #endif /* if WMOPS */
   }

    else {
      L40_var_out = L40_var1;

      for ( ; var2 > 0; var2--) {
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
         if( L40_var_out > 0x003fffffffff) {
#else
         if( L40_var_out > 0x003fffffffffLL) {
#endif
            L40_var_out = L40_OVERFLOW_OCCURED( L40_var_out);
            break;
         }

         else if ( L40_var_out < L40_constant) {
            L40_var_out = L40_UNDERFLOW_OCCURED( L40_var_out);
            break;
         }

         else {
            L40_var_out = L40_var_out << 1;
         }
      }
   }

   #if (WMOPS)
   multiCounter[currCounter].L40_set--;
   multiCounter[currCounter].L40_shl++;
   #endif /* if WMOPS */

   BASOP_CHECK(0, 0);

   return( L40_var_out);
}


/*****************************************************************************
 *
 *  Function Name : L40_shr
 *
 *  Purpose :
 *
 *    Arithmetically shifts right L40_var1 by var2 positions.
 *    - If var2 is positive, L40_var1 is shifted to the LSBits by (var2)
 *      positions with extension of the sign bit.
 *    - If var2 is negative, L40_var1 is shifted to the MSBits by (-var2)
 *      positions.
 *    Calls the macro L40_UNDERFLOW_OCCURED() in case of underflow on 40-bit.
 *    Calls the macro L40_OVERFLOW_OCCURED()  in case of overflow  on 40-bit.
 *
 *  Complexity weight : 1
 *
 *  Inputs :
 *
 *    L40_var1    40 bit long signed integer (Word40) whose value falls in the
 *                range : MIN_40 <= L40_var1 <= MAX_40.
 *
 *    var2        16 bit short signed integer (Word16) whose value falls in
 *                the range : MIN_16 <= var2 <= MAX_16.
 *
 *  Outputs :
 *
 *    none
 *
 *  Return Value :
 *
 *    L40_var_out 40 bit long signed integer (Word40) whose value falls in
 *                the range : MIN_40 <= L40_var_out <= MAX_40.
 *
 *****************************************************************************/
static __inline Word40 L40_shr( Word40 L40_var1, Word16 var2) {
   Word40 L40_var_out;

   if( var2 < 0) {
      var2 = -var2;
      L40_var_out  = L40_shl ( L40_var1, var2);

      #if (WMOPS)
      multiCounter[currCounter].L40_shl--;
      #endif /* if WMOPS */

   } else {
      L40_var_out = L40_var1 >> var2;

   }

   #if (WMOPS)
   multiCounter[currCounter].L40_shr++;
   #endif /* if WMOPS */

   return( L40_var_out);
}

/* end of file */
