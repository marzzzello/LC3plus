/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef __BASOP_UTIL_ROM_H__
#define __BASOP_UTIL_ROM_H__

#ifndef COUNT_ROM
#include "functions.h"
#endif
#include "typedef.h"

#define LD_INT_TAB_LEN 120
#define INV_TABLE_SIZE 256
#define SQRT_TABLE_SIZE 256

#ifndef CHEAP_NORM_SIZE
#define CHEAP_NORM_SIZE 161
#endif

#define MINSFTAB 7
#define MAXSFTAB 25

#define SHC(x) ((Word16)x)

/**
 * \brief  Lookup-Table for binary logarithm
 */
extern const Word16 ldCoeff[7];

/**
  \brief 	Lookup-Table for binary power algorithm
*/
extern const UWord32 exp2_tab_long[32];

/**
  \brief 	Lookup-Table for binary power algorithm
*/
extern const UWord32 exp2w_tab_long[32];

/**
  \brief 	Lookup-Table for binary power algorithm
*/
extern const UWord32 exp2x_tab_long[32];

/**
 * \brief 1/x, x=[0,1,2,3...]  table
 */
extern const Word16 InvIntTable[32];

/**
 * \ brief Sine tables
 */
extern const PWord16 SineTable480[241];
extern const PWord16 SineTable320[161];

/**
 * \ brief Lookup for sine tables and windows.
 */
void BASOP_getTables(const PWord16 **ptwiddle, const PWord16 **sin_twiddle, Word16 *sin_step, Word16 length);

extern const Word32 RealFFT20_twid[6];
extern const Word32 RealFFT32_twid[10];
extern const Word32 RealFFT40_twid[12];
extern const Word32 RealFFT60_twid[17];
extern const Word32 RealFFT64_twid[18];
extern const Word32 RealFFT80_twid[22];
extern const Word32 RealFFT96_twid[26];
extern const Word32 RealFFT128_twid[34];
extern const Word32 RealFFT192_twid[50];
extern const Word32 RealFFT256_twid[66];
extern const Word32 RealFFT384_twid[98];
extern const Word32 RealFFT512_twid[130];
extern const Word32 RealFFT768_twid[194];

extern const Word32 RotVector_32_32[2 * 20];
extern const Word32 RotVector_40_32[2 * 28];
extern const Word16 RotVector_320[2 * (320 - 20)];
extern const Word16 RotVector_360[2 * (360 - 30)];
extern const Word16 RotVector_480[2 * (480 - 30)];
extern const Word16 RotVector_32_8[2 * (256 - 32)];
extern const Word16 RotVector_32_12[2 * (384 - 32)];

extern const Word32 isqrt_table[128 + 2];

extern const Word32 Log2_16_table1[16];
extern const Word16 Log2_16_table2[16];

extern const Word32 InvLog2_16_table1[64];
extern const Word16 InvLog2_16_table2[64];

extern const UWord8  gf16_mult_table[256];
extern const UWord8  rs16_elp_deg2_table[256];
extern const UWord16 rs16_elp_deg3_table[256];

#endif
