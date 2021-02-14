/*
  ===========================================================================
   File: ENHUL32.H                                       v.0.6 - 11.Nov.2016
  ===========================================================================


  ============================================================================
*/


#ifndef _ENHUL32_H
#define _ENHUL32_H


/*****************************************************************************
 *
 *  Constants and Globals
 *
 *****************************************************************************/
#define ENHUL32               /* all  the  enhanced unsigned operators */
#define STL_TYPECASTS         /*  logical shift and bitwise manipulation functions              */
                              /*  algorithmically exact to existing signed L_lshr and L_lshr    */ 

/* #define DRAFT_UL_UPDATE  */  /* editorial corrections and speed improvements to UL_addNs and UL_subNs  */

#include "stl.h"


#ifndef UWord64 
#define UWord64 unsigned long long      /*  for local use inside UL_Mpy_32_*   */
#endif




#if (WMOPS)
#include "count.h"
extern BASIC_OP multiCounter[MAXCOUNTERS];  /* existing signed counters are reused for unsigedn operators */
extern int currCounter;
#endif /* if WMOPS */

/*****************************************************************************
 *
 *  Prototypes for enhanced unsigned 32 bit arithmetic operators
 *
 *****************************************************************************/
UWord32     UL_addNs(UWord32 a, UWord32 b, UWord16* wrap);  
#ifdef DRAFT_UL_UPDATE  
UWord32     UL_subNs(UWord32 a, UWord32 b, UWord16* sgn);  
#else
UWord32     UL_subNs(UWord32 a, UWord32 b, UWord16* wrap);  
#endif

UWord32  UL_Mpy_32_32(UWord32 a, UWord32 b);                            
void Mpy_32_32_uu( UWord32 a, UWord32 b, UWord32 *c_h, UWord32 *c_l);   /* does not saturate */
void Mpy_32_16_uu( UWord32 a, UWord16 b,UWord32 *c_h, UWord16 *c_l);    /* does not saturate   */

/*  Other  */
Word16  norm_ul (UWord32 UL_var1);  
UWord32 UL_deposit_l(UWord16);      /* deposit low without sign extension ) */


/*****************************************************************************
 *
 *  Inline Functions
 *
 *****************************************************************************/

#ifdef STL_TYPECASTS 
/*      (Reuse of existing signed STL "L" operators) with
        typecasting  to  make the resulting "UL" code a lot cleaner and more readable. */
UWord32 UL_lshl( UWord32 UL_var1, Word16 var2);
UWord32 UL_lshr( UWord32 UL_var1, Word16 var2);
UWord32 UL_and(UWord32 UL_var1, UWord32 UL_var2 );
UWord32 UL_or(UWord32 UL_var1, UWord32 UL_var2 );
UWord32 UL_xor(UWord32 UL_var1, UWord32 UL_var2 );
UWord32 UL_deposit_h(UWord16 uvar1);
UWord16 u_extract_h(UWord32 UL_var1);
UWord16 u_extract_l(UWord32 UL_var1);

/* enable convenient reuse of Non-saturating UL_subNs , UL_addNs  
   while "D"iscarding the sgn/wrap output flags */ 
UWord32 UL_subNsD(UWord32 UL_var1, UWord32 UL_var2 );
UWord32 UL_addNsD(UWord32 UL_var1, UWord32 UL_var2 );
#endif

#endif /*_ENHUL32_H*/

/* end of file */
