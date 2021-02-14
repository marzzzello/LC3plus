/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef __BASOP_UTIL_H__
#define __BASOP_UTIL_H__

#include "basop32.h"
#include "basop_mpy.h"
#include "defines.h"
#include "stl.h"
#include "typedef.h"
#include <assert.h>
#include <string.h>

#define _LONG long
#define _SHORT short
#ifdef _WIN32
#define _INT64 __int64
#else
#define _INT64 long long
#endif

#define WORD32_BITS 32
#define MAXVAL_WORD32 ((signed)0x7FFFFFFF)
#define MINVAL_WORD32 ((signed)0x80000000)
#define WORD32_FIX_SCALE ((_INT64)(1) << (WORD32_BITS - 1))

#define WORD16_BITS 16
#define MAXVAL_WORD16 (((signed)0x7FFFFFFF) >> 16)
#define MINVAL_WORD16 (((signed)0x80000000) >> 16)
#define WORD16_FIX_SCALE ((_INT64)(1) << (WORD16_BITS - 1))

/*!
  \def   Macro converts a Word32 fixed point to Word16 fixed point <1 with saturation
*/
#define WORD322WORD16(val)                                                                                             \
    ((((((val) >> (WORD32_BITS - WORD16_BITS - 1)) + 1) > (((_LONG)1 << WORD16_BITS) - 1)) && ((_LONG)(val) > 0))      \
         ? (Word16)(_SHORT)(((_LONG)1 << (WORD16_BITS - 1)) - 1)                                                       \
         : (Word16)(_SHORT)((((val) >> (WORD32_BITS - WORD16_BITS - 1)) + 1) >> 1))


/* Word16 Packed Type */
typedef struct
{
    struct
    {
        Word16 re;
        Word16 im;
    } v;
} PWord16;

#define cast16 move16

#define LD_DATA_SCALE (6)
#define LD_DATA_SHIFT_I5 (7)

/************************************************************************/
/*!
  \brief   Calculate the squareroot of a number given by mantissa and exponent

  Mantissa is in 16/32-bit-fractional format with values between 0 and 1. <br>
  For *norm versions mantissa has to be between 0.5 and 1. <br>
  The base for the exponent is 2.  Example:  \f$  a = a\_m * 2^{a\_e}  \f$<br>
  The exponent is addressed via pointers and will be overwritten with the result.
*/

Word16 Sqrt16(                  /*!< output mantissa */
              Word16  mantissa, /*!< input mantissa */
              Word16 *exponent  /*!< pointer to exponent */
);

Word16 ISqrt16(                  /*!< output mantissa */
               Word16  mantissa, /*!< input mantissa */
               Word16 *exponent  /*!< pointer to exponent */
);

/*****************************************************************************/
/*!
  \brief   Calculate the inverse of a number given by mantissa and exponent

  Mantissa is in 16-bit-fractional format with values between 0 and 1. <br>
  The base for the exponent is 2.  Example:  \f$  a = a\_m * 2^{a\_e}  \f$<br>
  The operand is addressed via pointers and will be overwritten with the result.

  The function uses a table lookup and a newton iteration.
*/
Word16 Inv16(                  /*!< output mantissa */
             Word16  mantissa, /*!< input mantissa */
             Word16 *exponent  /*!< pointer to exponent */
);

/********************************************************************/
/*!
  \brief   Calculates the scalefactor needed to normalize input array

    The scalefactor needed to normalize the Word16 input array is returned <br>
    If the input array contains only '0', a scalefactor 0 is returned <br>
    Scaling factor is determined wrt a normalized target x: 16384 <= x <= 32767 for positive x <br>
    and   -32768 <= x <= -16384 for negative x
*/

Word16 getScaleFactor16(                     /* o: measured headroom in range [0..15], 0 if all x[i] == 0 */
                        const Word16 *x,     /* i: array containing 16-bit data */
                        const Word16  len_x); /* i: length of the array to scan  */

/********************************************************************/
/*!
  \brief   Calculates the scalefactor needed to normalize input array

    The scalefactor needed to normalize the Word16 input array is returned <br>
    If the input array contains only '0', a scalefactor 16 is returned <br>
    Scaling factor is determined wrt a normalized target x: 16384 <= x <= 32767 for positive x <br>
    and   -32768 <= x <= -16384 for negative x
*/

Word16 getScaleFactor16_0(                     /* o: measured headroom in range [0..15], 16 if all x[i] == 0 */
                          const Word16 *x,     /* i: array containing 16-bit data */
                          const Word16  len_x); /* i: length of the array to scan  */

/********************************************************************/
/*!
  \brief   Calculates the scalefactor needed to normalize input array

    The scalefactor needed to normalize the Word32 input array is returned <br>
    If the input array contains only '0', a scalefactor 0 is returned <br>
    Scaling factor is determined wrt a normalized target x: 1073741824 <= x <= 2147483647 for positive x <br>
    and   -2147483648 <= x <= -1073741824 for negative x
*/

Word16 getScaleFactor32(                     /* o: measured headroom in range [0..31], 0 if all x[i] == 0 */
                        const Word32 *x,     /* i: array containing 32-bit data */
                        const Word16  len_x); /* i: length of the array to scan  */

/********************************************************************/
/*!
  \brief   Calculates the scalefactor needed to normalize input array

    The scalefactor needed to normalize the Word32 input array is returned <br>
    If the input array contains only '0', a scalefactor 32 is returned <br>
    Scaling factor is determined wrt a normalized target x: 1073741824 <= x <= 2147483647 for positive x <br>
    and   -2147483648 <= x <= -1073741824 for negative x
*/

Word16 getScaleFactor32_0(                     /* o: measured headroom in range [0..31], 32 if all x[i] == 0 */
                          const Word32 *x,     /* i: array containing 32-bit data */
                          const Word16  len_x); /* i: length of the array to scan  */

/****************************************************************************/
/*!
  \brief   Does fractional integer division of Word32 arg1 by Word16 arg2


  \return fractional Word16 integer z = arg1(32bits)/arg2(16bits) , z not normalized
*/
Word16 BASOP_Util_Divide3216_Scale(Word32  x,  /*!< i  : Numerator  */
                                   Word16  y,  /*!< i  : Denominator*/
                                   Word16 *s); /*!< o  : Additional scalefactor difference*/

/****************************************************************************/
/*!
  \brief   Does fractional division of Word16 arg1 by Word16 arg2


  \return fractional Q15 Word16 z = arg1(Q15)/arg2(Q15)  with scaling s
*/
Word16 BASOP_Util_Divide1616_Scale(Word16  x,  /*!< i  : Numerator*/
                                   Word16  y,  /*!< i  : Denominator*/
                                   Word16 *s); /*!< o  : Additional scalefactor difference*/

/************************************************************************/
/*!
  \brief 	Binary logarithm with 7 iterations

  \param   x

  \return log2(x)/64
 */
/************************************************************************/
Word32 BASOP_Util_Log2(Word32 x);

/************************************************************************/
/*!
  \brief 	Binary power

  Date: 06-JULY-2012 Arthur Tritthart, IIS Fraunhofer Erlangen

  Version with 3 table lookup and 1 linear interpolations

  Algorithm: compute power of 2, argument x is in Q7.25 format
             result = 2^(x/64)
             We split exponent (x/64) into 5 components:
             integer part:      represented by b31..b25  (exp)
             fractional part 1: represented by b24..b20  (lookup1)
             fractional part 2: represented by b19..b15  (lookup2)
             fractional part 3: represented by b14..b10  (lookup3)
             fractional part 4: represented by b09..b00  (frac)
             => result = (lookup1*lookup2*(lookup3+C1*frac)<<3)>>exp

  Due to the fact, that all lookup values contain a factor 0.5
  the result has to be shifted by 3 to the right also.
  Table exp2_tab_long contains the log2 for 0 to 1.0 in steps
  of 1/32, table exp2w_tab_long the log2 for 0 to 1/32 in steps
  of 1/1024, table exp2x_tab_long the log2 for 0 to 1/1024 in
  steps of 1/32768. Since the 2-logarithm of very very small
  negative value is rather linear, we can use interpolation.

  Limitations:

  For x <= 0, the result is fractional positive
  For x > 0, the result is integer in range 1...7FFF.FFFF
  For x < -31/64, we have to clear the result
  For x = 0, the result is ~1.0 (0x7FFF.FFFF)
  For x >= 31/64, the result is 0x7FFF.FFFF

  \param  x

  \return pow(2,(x/64))
 */
/************************************************************************/
Word32 BASOP_Util_InvLog2(Word32 x);

/**
 * \brief Compute dot product of 1 32 bit vectors with itself
 * \param x input vector 1
 * \param headroom amount of headroom bits the input vector
 * \param length the length of the input vector
 * \param result_e pointer to where the exponent of the result will be stored into
 * \return the dot product of x and x.
 */
Word32 Norm32Norm(const Word32 *x, const Word16 headroom, const Word16 length, Word16 *result_e);

/*!**********************************************************************
   \brief   Add two values given by mantissa and exponent.

   Mantissas are in 32-bit-fractional format with values between 0 and 1. <br>
   The base for exponents is 2.  Example:  \f$  a = a\_m * 2^{a\_e}  \f$<br>

************************************************************************/
Word32 BASOP_Util_Add_Mant32Exp /*!< o: normalized result mantissa */
    (Word32  a_m,               /*!< i: Mantissa of 1st operand a  */
     Word16  a_e,               /*!< i: Exponent of 1st operand a  */
     Word32  b_m,               /*!< i: Mantissa of 2nd operand b  */
     Word16  b_e,               /*!< i: Exponent of 2nd operand b  */
     Word16 *ptr_e);            /*!< o: exponent of result         */
                                /*!**********************************************************************
                                   \brief   Returns the comparison result of two normalized values given by mantissa and exponent.
                                            return value: -1: a < b, 0: a == b, 1; a > b

                                   Mantissas are in 32-bit-fractional format with values between 0 and 1. <br>
                                   The base for exponents is 2.  Example:  \f$  a = a\_m * 2^{a\_e}  \f$<br>

                                ************************************************************************/
Word16 BASOP_Util_Cmp_Mant32Exp /*!< o: flag: result of comparison */
    (Word32 a_m,                /*!< i: Mantissa of 1st operand a  */
     Word16 a_e,                /*!< i: Exponent of 1st operand a  */
     Word32 b_m,                /*!< i: Mantissa of 2nd operand b  */
     Word16 b_e);               /*!< i: Exponent of 2nd operand b  */

/* compare two positive normalized 16 bit mantissa/exponent values */
/* return value: positive if first value greater, negative if second value greater, zero if equal */
Word16 compMantExp16Unorm(Word16 m1, Word16 e1, Word16 m2, Word16 e2);

void Scale_sig(Word16       x[], /* i/o: signal to scale                 Qx        */
               const Word16 lg,  /* i  : size of x[]                     Q0        */
               const Word16 exp0 /* i  : exponent: x = round(x << exp)   Qx ?exp  */
);
void Copy_Scale_sig(const Word16 x[], /* i  : signal to scale input           Qx        */
                    Word16       y[], /* o  : scaled signal output            Qx        */
                    const Word16 lg,  /* i  : size of x[]                     Q0        */
                    const Word16 exp0 /* i  : exponent: x = round(x << exp)   Qx ?exp  */
);

Word32 Isqrt(Word32  x,  /* (i)   Q31: normalized value (1.0 > x >= 0.5) */
             Word16 *x_e /* (i/o) Q0 : pointer to exponent */
);

Word16 BASOP_Util_Log2_16(Word32 x, Word16 x_e);

Word16 BASOP_Util_InvLog2_16(Word16 x, Word16 *y_e);

#ifdef USE_KISS_FFT
void fft16(Word32 *re, Word32 *im, Word16 s);
#endif

#define BASOP_CFFT_MAX_LENGTH 384
void BASOP_cfft(Word32 *re, Word32 *im, Word16 sizeOfFft, Word16 s, Word16 *scale, Word32 *x);
void BASOP_rfftN(Word32 *re, Word16 sizeOfFft, Word16 *scale, Word8 *scratchBuffer);
void BASOP_irfftN(Word32 *re, Word16 sizeOfFft, Word16 *scale, Word8 *scratchBuffer);

#if WMOPS
extern BASIC_OP multiCounter[MAXCOUNTERS];
extern int      currCounter;
#endif

static __inline void basop_memcpy(void *dst, const void *src, size_t n)
{
#if WMOPS
    multiCounter[currCounter].move16 += (UWord32)n / 2;
#endif
    /* check for overlapping memory */
    assert((const char *)src + n <= (char *)dst || (char *)dst + n <= (const char *)src);
    memcpy(dst, src, n);
}

static __inline void basop_memmove(void *dst, const void *src, size_t n)
{
#if WMOPS
    multiCounter[currCounter].move16 += (UWord32)n / 2;
#endif
    memmove(dst, src, n);
}

static __inline void basop_memset(void *dst, int val, size_t n)
{
#if WMOPS
    multiCounter[currCounter].move16 += (UWord32)n / 2;
#endif
    memset(dst, val, n);
}

/* Macros around Dyn_Mem that don't require duplicate declarations. */
#ifdef DYNMEM_COUNT
/* older visual studio doesn't have __func__ */
#if defined _MSC_VER && _MSC_VER < 1900
#define __func__ __FUNCTION__
#endif
#define Dyn_Mem_Deluxe_In(...) __VA_ARGS__ Dyn_Mem_In(__func__, sizeof(struct {__VA_ARGS__}))
#define Dyn_Mem_Deluxe_Out() Dyn_Mem_Out()
#else
#define Dyn_Mem_Deluxe_In(...) __VA_ARGS__
#define Dyn_Mem_Deluxe_Out()
#endif

#endif /* __BASOP_UTIL_H__ */
