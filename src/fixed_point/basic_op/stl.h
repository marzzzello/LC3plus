/*
  ===========================================================================
   File: STL.H                                           v.2.3 - 30.Nov.2009
  ===========================================================================

            ITU-T STL  BASIC OPERATORS

            MAIN      HEADER      FILE

   History:
   07 Nov 04   v2.0     Incorporation of new 32-bit / 40-bit / control
                        operators for the ITU-T Standard Tool Library as 
                        described in Geneva, 20-30 January 2004 WP 3/16 Q10/16
                        TD 11 document and subsequent discussions on the
                        wp3audio@yahoogroups.com email reflector.
   March 06   v2.1      Changed to improve portability.
   31 Mar 15  v2.1E     Clarified usage of operators needed for the EVS codec.

  ============================================================================
*/


#ifndef _STL_H
#define _STL_H

/* This is Defined right here, it should be put in an include file
 * that is included from every file. So far, only "stl.h" fit this
 * requirement as "options.h" and "options.h" are not included by
 * all code files. Also, these definitions are not enclosed by
 * #if DEBUG because in some files, "stl.h" is included before
 * DEBUG had been defined. This would lead to these macros not
 * being defined and the Stack Counting would be disabled for some
 * files. The Stack Counting would still work but some functions
 * would be missing from the tree. */
extern int check_stack(const char*, const char*);
#define STACK_DEPTH_FCT_CALL   (BASOP_push_wmops(__FUNCTION__), check_stack(   __FILE__, __FUNCTION__))
#define STACK_DEPTH_FCT_RETURN (BASOP_pop_wmops(), check_stack("-"__FILE__, __FUNCTION__))

#include <stdlib.h> /* for size_t */
extern void* check_alloc_in(const char*, const char *, size_t, size_t);
extern void  check_alloc_out(void*);
extern void  check_alloc_exit(void);

#define MALLOC_FCT_CALL(n1)     check_alloc_in("m"__FILE__, __FUNCTION__, n1, 0)
#define CALLOC_FCT_CALL(n1, n2) check_alloc_in("c"__FILE__, __FUNCTION__, n1, n2)
#define   FREE_FCT_CALL(ptr)    check_alloc_out(ptr)

/* both ALLOW_40bits and ALLOW_ENH_UL32 shall be enabled for the EVS codec.  */
#define ALLOW_40bits         /* allow 32x16 and 32x32 multiplications */
#define ALLOW_ENH_UL32       /* allow enhanced  unsigned 32bit  operators */


#include "typedef.h"
#include "basop32.h" 
#include "count.h"
#include "move.h"
#include "control.h"
#include "enh1632.h" 

#if defined (ALLOW_40bits)
#include "enh40.h"
#endif

#if defined (ALLOW_ENH_UL32)
#include "enhUL32.h"
#endif


#endif /* ifndef _STL_H */


/* end of file */
