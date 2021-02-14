/*
  ===========================================================================
   File: ENH40.H                                         v.2.3 - 30.Nov.2009
  ===========================================================================

            ITU-T  STL  BASIC OPERATORS

            40-BIT ARITHMETIC OPERATORS

   History:
   07 Nov 04   v2.0     Incorporation of new 32-bit / 40-bit / control
                        operators for the ITU-T Standard Tool Library as 
                        described in Geneva, 20-30 January 2004 WP 3/16 Q10/16
                        TD 11 document and subsequent discussions on the
                        wp3audio@yahoogroups.com email reflector.

   March 06    v2.1     Changed to improve portability.

   31 Mar 15   v2.1E    Removal of operators not used in the EVS codec.

  ============================================================================
*/


#ifndef _ENH40_H
#define _ENH40_H

#include "stl.h"
 
 /*****************************************************************************
 *
 *  Prototypes for enhanced 40 bit arithmetic operators
 *
 *****************************************************************************/

void Mpy_32_16_ss( Word32 L_var1, Word16 var2,   Word32 *L_varout_h, UWord16 *varout_l);
void Mpy_32_32_ss( Word32 L_var1, Word32 L_var2, Word32 *L_varout_h, UWord32 *L_varout_l);

#endif /*_ENH40_H*/


/* end of file */


