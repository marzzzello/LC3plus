/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"


#ifndef USE_KISS_FFT
#define SCALEFACTORN2 (3)
#define SCALEFACTOR4 (3)
#define SCALEFACTOR5 (4)
#define SCALEFACTOR6 (4)
#define SCALEFACTOR8 (4)
#define SCALEFACTOR15 (5)
#define SCALEFACTOR30_1 (5)
#define SCALEFACTOR30_2 (1)
#define SCALEFACTOR32_1 (5)
#define SCALEFACTOR32_2 (1)
#define Mpy_32_xx Mpy_32_16
#endif

#define SCALEFACTOR10 5
#define SCALEFACTOR16 5
#define SCALEFACTOR20 5
#define SCALEFACTOR30 6
#define SCALEFACTOR32 6
#define SCALEFACTOR40 7
#define SCALEFACTOR48 8
#define SCALEFACTOR60 7
#define SCALEFACTOR64 7
#define SCALEFACTOR80 8
#define SCALEFACTOR90 9
#define SCALEFACTOR96 9
#define SCALEFACTOR120 8
#define SCALEFACTOR128 8
#define SCALEFACTOR160 8
#define SCALEFACTOR180 10
#define SCALEFACTOR192 10
#define SCALEFACTOR240 9
#define SCALEFACTOR256 9
#define SCALEFACTOR384 11

#ifndef USE_KISS_FFT

#define FFTC(x) WORD322WORD16((Word32)x)

#define C31 (FFTC(0x91261468)) /* FL2WORD32( -0.86602540) -sqrt(3)/2 */

#define C51 (FFTC(0x79bc3854)) /* FL2WORD32( 0.95105652)   */
#define C52 (FFTC(0x9d839db0)) /* FL2WORD32(-1.53884180/2) */
#define C53 (FFTC(0xd18053ce)) /* FL2WORD32(-0.36327126)   */
#define C54 (FFTC(0x478dde64)) /* FL2WORD32( 0.55901699)   */
#define C55 (FFTC(0xb0000001)) /* FL2WORD32(-1.25/2)       */

#define C81 (FFTC(0x5a82799a)) /* FL2WORD32( 7.071067811865475e-1) */
#define C82 (FFTC(0xa57d8666)) /* FL2WORD32(-7.071067811865475e-1) */

#define C161 (FFTC(0x5a82799a)) /* FL2WORD32( 7.071067811865475e-1)  INV_SQRT2    */
#define C162 (FFTC(0xa57d8666)) /* FL2WORD32(-7.071067811865475e-1) -INV_SQRT2    */

#define C163 (FFTC(0x7641af3d)) /* FL2WORD32( 9.238795325112867e-1)  COS_PI_DIV8  */
#define C164 (FFTC(0x89be50c3)) /* FL2WORD32(-9.238795325112867e-1) -COS_PI_DIV8  */

#define C165 (FFTC(0x30fbc54d)) /* FL2WORD32( 3.826834323650898e-1)  COS_3PI_DIV8 */
#define C166 (FFTC(0xcf043ab3)) /* FL2WORD32(-3.826834323650898e-1) -COS_3PI_DIV8 */

#define C51_32 (0x79bc3854) /* FL2WORD32( 0.95105652)   */
#define C52_32 (0x9d839db0) /* FL2WORD32(-1.53884180/2) */
#define C53_32 (0xd18053ce) /* FL2WORD32(-0.36327126)   */
#define C54_32 (0x478dde64) /* FL2WORD32( 0.55901699)   */
#define C55_32 (0xb0000001) /* FL2WORD32(-1.25/2)       */

#define C61_32 (0x6ed9eba1)

#define C81_32 (0x5a82799a) /* FL2WORD32( 7.071067811865475e-1) */
#define C82_32 (0xa57d8666) /* FL2WORD32(-7.071067811865475e-1) */

#define Mpy3_0(s12, s13, s14, s15, t0, t1, t2, t3)                                                                     \
    do                                                                                                                 \
    {                                                                                                                  \
        s12 = Mpy_32_32(L_add(t0, t2), C81_32);                                                                        \
        s14 = Mpy_32_32(L_sub(t0, t2), C81_32);                                                                        \
        s13 = Mpy_32_32(L_sub(t3, t1), C81_32);                                                                        \
        s15 = Mpy_32_32(L_add(t1, t3), C82_32);                                                                        \
    } while (0)

#define cplxMpy3_0(a, b, c, d)                                                                                         \
    do                                                                                                                 \
    {                                                                                                                  \
        as = L_shr_pos(a, 1);                                                                                          \
        bs = L_shr_pos(b, 1);                                                                                          \
        a  = L_sub(Mpy_32_32(as, c), Mpy_32_32(bs, d));                                                                \
        b  = L_add(Mpy_32_32(as, d), Mpy_32_32(bs, c));                                                                \
    } while (0)

#define cplxMpy4_4_0(re, im, a, b, c, d)                                                                               \
    re = L_shr(L_sub(Mpy_32_xx(a, c), Mpy_32_xx(b, d)), SCALEFACTOR60 - SCALEFACTOR15);                                \
    im = L_shr(L_add(Mpy_32_xx(a, d), Mpy_32_xx(b, c)), SCALEFACTOR60 - SCALEFACTOR15);

#define cplxMpy4_4_1(re, im, a, b)                                                                                     \
    re = L_shr(a, SCALEFACTOR60 - SCALEFACTOR15);                                                                      \
    im = L_shr(b, SCALEFACTOR60 - SCALEFACTOR15);


#define cplxMpy4_8_0(re, im, a, b, c, d)                                                                               \
    do                                                                                                                 \
    {                                                                                                                  \
        re = L_shr_pos(L_sub(Mpy_32_xx(a, c), Mpy_32_xx(b, d)), 1);                                                    \
        im = L_shr_pos(L_add(Mpy_32_xx(a, d), Mpy_32_xx(b, c)), 1);                                                    \
    } while (0)


#define cplxMpy4_8_1(re, im, a, b)                                                                                     \
    do                                                                                                                 \
    {                                                                                                                  \
        re = L_shr_pos(a, 1);                                                                                          \
        im = L_shr_pos(b, 1);                                                                                          \
    } while (0)


#define cplxMpy4_8_2(re, im, a, b, c, d)                                                                               \
    do                                                                                                                 \
    {                                                                                                                  \
        re = L_shr_pos(L_add(Mpy_32_32(a, c), Mpy_32_32(b, d)), 1);                                                    \
        im = L_shr_pos(L_sub(Mpy_32_32(b, c), Mpy_32_32(a, d)), 1);                                                    \
    } while (0)


#define cplxMpy4_12_0(re, im, a, b, c, d)                                                                              \
    do                                                                                                                 \
    {                                                                                                                  \
        re = L_sub(Mpy_32_xx(a, c), Mpy_32_xx(b, d));                                                                  \
    move32();                                                                                                     \
        im = L_add(Mpy_32_xx(a, d), Mpy_32_xx(b, c));                                                                  \
        move32();                                                                                                     \
    } while (0)

#define cplxMpy4_12_1(re, im, a, b)                                                                                    \
    do                                                                                                                 \
    {                                                                                                                  \
        re = a;                                                                                                        \
    move32();                                                                                                     \
        im = b;                                                                                                        \
        move32();                                                                                                     \
    } while (0)

static void fft4(Word32 *x)
{
    Dyn_Mem_Deluxe_In(
        Word32 x0, x1, x2, x3, x4, x5, x6, x7;
        Word32 t0, t1, t2, t3, t4, t5, t6, t7;
    );

    x0 = L_shr_pos(x[0], SCALEFACTOR4);
    x1 = L_shr_pos(x[1], SCALEFACTOR4);
    x2 = L_shr_pos(x[2], SCALEFACTOR4);
    x3 = L_shr_pos(x[3], SCALEFACTOR4);
    x4 = L_shr_pos(x[4], SCALEFACTOR4);
    x5 = L_shr_pos(x[5], SCALEFACTOR4);
    x6 = L_shr_pos(x[6], SCALEFACTOR4);
    x7 = L_shr_pos(x[7], SCALEFACTOR4);

    /* Pre-additions */
    t0 = L_add(x0, x4);
    t2 = L_sub(x0, x4);
    t1 = L_add(x1, x5);
    t3 = L_sub(x1, x5);
    t4 = L_add(x2, x6);
    t7 = L_sub(x2, x6);
    t5 = L_add(x7, x3);
    t6 = L_sub(x7, x3);

    /* Post-additions */
    x[0] = L_add(t0, t4);
    x[1] = L_add(t1, t5);
    x[2] = L_sub(t2, t6);
    x[3] = L_sub(t3, t7);
    x[4] = L_sub(t0, t4);
    x[5] = L_sub(t1, t5);
    x[6] = L_add(t2, t6);
    x[7] = L_add(t3, t7);

    Dyn_Mem_Deluxe_Out();
}

/**
 * \brief    Function performs a complex 5-point FFT
 *           The FFT is performed inplace. The result of the FFT
 *           is scaled by SCALEFACTOR5 bits.
 *
 *
 * \param    [i/o] re    real input / output
 * \param    [i/o] im    imag input / output
 * \param    [i  ] s     stride real and imag input / output
 *
 * \return   void
 */


static void fft5(Word32 *re, Word32 *im, Word16 s)
{
    Dyn_Mem_Deluxe_In(
        Word32 x0, x1, x2, x3, x4;
        Word32 r1, r2, r3, r4;
        Word32 s1, s2, s3, s4;
        Word32 t;
    );

    /* real part */
    x0 = L_shr_pos(re[s * 0], SCALEFACTOR5);
    x1 = L_shr_pos(re[s * 1], SCALEFACTOR5);
    x2 = L_shr_pos(re[s * 2], SCALEFACTOR5);
    x3 = L_shr_pos(re[s * 3], SCALEFACTOR5);
    x4 = L_shr_pos(re[s * 4], SCALEFACTOR5);

    r1    = L_add(x1, x4);
    r4    = L_sub(x1, x4);
    r3    = L_add(x2, x3);
    r2    = L_sub(x2, x3);
    t     = Mpy_32_32(L_sub(r1, r3), C54_32);
    r1    = L_add(r1, r3);
    re[0] = L_add(x0, r1); move32();
    /* Bit shift left because of the constant C55 which was scaled with the factor 0.5 because of the representation of
       the values as fracts */
    r1 = L_add(re[0], (L_shl_pos(Mpy_32_32(r1, C55_32), 1)));
    r3 = L_sub(r1, t);
    r1 = L_add(r1, t);
    t  = Mpy_32_32(L_add(r4, r2), C51_32);
    /* Bit shift left because of the constant C55 which was scaled with the factor 0.5 because of the representation of
       the values as fracts */
    r4 = L_add(t, L_shl_pos(Mpy_32_32(r4, C52_32), 1));
    r2 = L_add(t, Mpy_32_32(r2, C53_32));

    /* imaginary part */
    x0 = L_shr_pos(im[s * 0], SCALEFACTOR5);
    x1 = L_shr_pos(im[s * 1], SCALEFACTOR5);
    x2 = L_shr_pos(im[s * 2], SCALEFACTOR5);
    x3 = L_shr_pos(im[s * 3], SCALEFACTOR5);
    x4 = L_shr_pos(im[s * 4], SCALEFACTOR5);

    s1    = L_add(x1, x4);
    s4    = L_sub(x1, x4);
    s3    = L_add(x2, x3);
    s2    = L_sub(x2, x3);
    t     = Mpy_32_32(L_sub(s1, s3), C54_32);
    s1    = L_add(s1, s3);
    im[0] = L_add(x0, s1); move32();
    /* Bit shift left because of the constant C55 which was scaled with the factor 0.5 because of the representation of
       the values as fracts */
    s1 = L_add(im[0], L_shl_pos(Mpy_32_32(s1, C55_32), 1));
    s3 = L_sub(s1, t);
    s1 = L_add(s1, t);
    t  = Mpy_32_32(L_add(s4, s2), C51_32);
    /* Bit shift left because of the constant C55 which was scaled with the factor 0.5 because of the representation of
       the values as fracts */
    s4 = L_add(t, L_shl_pos(Mpy_32_32(s4, C52_32), 1));
    s2 = L_add(t, Mpy_32_32(s2, C53_32));

    /* combination */
    re[s * 1] = L_add(r1, s2); move32();
    re[s * 4] = L_sub(r1, s2); move32();
    re[s * 2] = L_sub(r3, s4); move32();
    re[s * 3] = L_add(r3, s4); move32();

    im[s * 1] = L_sub(s1, r2); move32();
    im[s * 4] = L_add(s1, r2); move32();
    im[s * 2] = L_add(s3, r4); move32();
    im[s * 3] = L_sub(s3, r4); move32();

    Dyn_Mem_Deluxe_Out();
}


/**
 * \brief    Function performs a complex 6-point FFT
 *           The FFT is performed inplace. The result of the FFT
 *           is scaled by SCALEFACTOR6 bits.
 *
 *
 * \param    [i/o] re    real input / output
 * \param    [i/o] im    imag input / output
 * \param    [i  ] st     stride real and imag input / output
 *
 * \return   void
 */



/**
 * \brief    Function performs a complex 8-point FFT
 *           The FFT is performed inplace. The result of the FFT
 *           is scaled by SCALEFACTOR8 bits.
 *
 *           WOPS with 32x16 bit multiplications: 108 cycles
 *
 * \param    [i/o] re    real input / output
 * \param    [i/o] im    imag input / output
 * \param    [i  ] s     stride real and imag input / output
 *
 * \return   void
 */


static void fft8(Word32 *re, Word32 *im, Word16 s)
{
    Dyn_Mem_Deluxe_In(
        Word32 x00, x01, x02, x03, x04, x05, x06, x07;
        Word32 x08, x09, x10, x11, x12, x13, x14, x15;
        Word32 t00, t01, t02, t03, t04, t05, t06, t07;
        Word32 t08, t09, t10, t11, t12, t13, t14, t15;
        Word32 s00, s01, s02, s03, s04, s05, s06, s07;
        Word32 s08, s09, s10, s11, s12, s13, s14, s15;
    );

    /* Pre-additions */

    x00 = L_shr_pos(re[s * 0], SCALEFACTOR8);
    x01 = L_shr_pos(im[s * 0], SCALEFACTOR8);
    x02 = L_shr_pos(re[s * 1], SCALEFACTOR8);
    x03 = L_shr_pos(im[s * 1], SCALEFACTOR8);
    x04 = L_shr_pos(re[s * 2], SCALEFACTOR8);
    x05 = L_shr_pos(im[s * 2], SCALEFACTOR8);
    x06 = L_shr_pos(re[s * 3], SCALEFACTOR8);
    x07 = L_shr_pos(im[s * 3], SCALEFACTOR8);
    x08 = L_shr_pos(re[s * 4], SCALEFACTOR8);
    x09 = L_shr_pos(im[s * 4], SCALEFACTOR8);
    x10 = L_shr_pos(re[s * 5], SCALEFACTOR8);
    x11 = L_shr_pos(im[s * 5], SCALEFACTOR8);
    x12 = L_shr_pos(re[s * 6], SCALEFACTOR8);
    x13 = L_shr_pos(im[s * 6], SCALEFACTOR8);
    x14 = L_shr_pos(re[s * 7], SCALEFACTOR8);
    x15 = L_shr_pos(im[s * 7], SCALEFACTOR8);

    t00 = L_add(x00, x08);
    t02 = L_sub(x00, x08);
    t01 = L_add(x01, x09);
    t03 = L_sub(x01, x09);
    t04 = L_add(x02, x10);
    t06 = L_sub(x02, x10);
    t05 = L_add(x03, x11);
    t07 = L_sub(x03, x11);
    t08 = L_add(x04, x12);
    t10 = L_sub(x04, x12);
    t09 = L_add(x05, x13);
    t11 = L_sub(x05, x13);
    t12 = L_add(x06, x14);
    t14 = L_sub(x06, x14);
    t13 = L_add(x07, x15);
    t15 = L_sub(x07, x15);

    /* Pre-additions and core multiplications */

    s00 = L_add(t00, t08);
    s04 = L_sub(t00, t08);
    s01 = L_add(t01, t09);
    s05 = L_sub(t01, t09);
    s08 = L_sub(t02, t11);
    s10 = L_add(t02, t11);
    s09 = L_add(t03, t10);
    s11 = L_sub(t03, t10);
    s02 = L_add(t04, t12);
    s07 = L_sub(t04, t12);
    s03 = L_add(t05, t13);
    s06 = L_sub(t13, t05);

    t01 = L_add(t06, t14);
    t02 = L_sub(t06, t14);
    t00 = L_add(t07, t15);
    t03 = L_sub(t07, t15);

    s12 = Mpy_32_xx(L_add(t00, t02), C81);
    s14 = Mpy_32_xx(L_sub(t00, t02), C81);
    s13 = Mpy_32_xx(L_sub(t03, t01), C81);
    s15 = Mpy_32_xx(L_add(t01, t03), C82);

    /* Post-additions */

    re[s * 0] = L_add(s00, s02); move32();
    re[s * 4] = L_sub(s00, s02); move32();
    im[s * 0] = L_add(s01, s03); move32();
    im[s * 4] = L_sub(s01, s03); move32();
    re[s * 2] = L_sub(s04, s06); move32();
    re[s * 6] = L_add(s04, s06); move32();
    im[s * 2] = L_sub(s05, s07); move32();
    im[s * 6] = L_add(s05, s07); move32();
    re[s * 3] = L_add(s08, s14); move32();
    re[s * 7] = L_sub(s08, s14); move32();
    im[s * 3] = L_add(s09, s15); move32();
    im[s * 7] = L_sub(s09, s15); move32();
    re[s * 1] = L_add(s10, s12); move32();
    re[s * 5] = L_sub(s10, s12); move32();
    im[s * 1] = L_add(s11, s13); move32();
    im[s * 5] = L_sub(s11, s13); move32();

    Dyn_Mem_Deluxe_Out();
}

/**
 * \brief    Function performs a complex 10-point FFT
 *           The FFT is performed inplace. The result of the FFT
 *           is scaled by SCALEFACTOR10 bits.
 *
 *           WOPS with 32x16 bit multiplications:  196 cycles
 *
 * \param    [i/o] re    real input / output
 * \param    [i/o] im    imag input / output
 * \param    [i  ] s     stride real and imag input / output
 *
 * \return   void
 */


static void fft10(Word32 *re, Word32 *im, Word16 s)
{
    Dyn_Mem_Deluxe_In(
        Word32 t;
        Word32 x0, x1, x2, x3, x4;
        Word32 r1, r2, r3, r4;
        Word32 s1, s2, s3, s4;
        Word32 y00, y01, y02, y03, y04, y05, y06, y07, y08, y09;
        Word32 y10, y11, y12, y13, y14, y15, y16, y17, y18, y19;
    );

    /* 2 fft5 stages */

    /* real part */
    x0 = L_shr_pos(re[s * 0], SCALEFACTOR10);
    x1 = L_shr_pos(re[s * 2], SCALEFACTOR10);
    x2 = L_shr_pos(re[s * 4], SCALEFACTOR10);
    x3 = L_shr_pos(re[s * 6], SCALEFACTOR10);
    x4 = L_shr_pos(re[s * 8], SCALEFACTOR10);

    r1  = L_add(x3, x2);
    r4  = L_sub(x3, x2);
    r3  = L_add(x1, x4);
    r2  = L_sub(x1, x4);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y00 = L_add(x0, r1);
    r1  = L_add(y00, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    x0 = L_shr_pos(im[s * 0], SCALEFACTOR10);
    x1 = L_shr_pos(im[s * 2], SCALEFACTOR10);
    x2 = L_shr_pos(im[s * 4], SCALEFACTOR10);
    x3 = L_shr_pos(im[s * 6], SCALEFACTOR10);
    x4 = L_shr_pos(im[s * 8], SCALEFACTOR10);

    s1  = L_add(x3, x2);
    s4  = L_sub(x3, x2);
    s3  = L_add(x1, x4);
    s2  = L_sub(x1, x4);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y01 = L_add(x0, s1);
    s1  = L_add(y01, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y04 = L_add(r1, s2);
    y16 = L_sub(r1, s2);
    y08 = L_sub(r3, s4);
    y12 = L_add(r3, s4);

    y05 = L_sub(s1, r2);
    y17 = L_add(s1, r2);
    y09 = L_add(s3, r4);
    y13 = L_sub(s3, r4);

    /* real part */
    x0 = L_shr_pos(re[s * 5], SCALEFACTOR10);
    x1 = L_shr_pos(re[s * 1], SCALEFACTOR10);
    x2 = L_shr_pos(re[s * 3], SCALEFACTOR10);
    x3 = L_shr_pos(re[s * 7], SCALEFACTOR10);
    x4 = L_shr_pos(re[s * 9], SCALEFACTOR10);

    r1  = L_add(x1, x4);
    r4  = L_sub(x1, x4);
    r3  = L_add(x3, x2);
    r2  = L_sub(x3, x2);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y02 = L_add(x0, r1);
    r1  = L_add(y02, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    x0 = L_shr_pos(im[s * 5], SCALEFACTOR10);
    x1 = L_shr_pos(im[s * 1], SCALEFACTOR10);
    x2 = L_shr_pos(im[s * 3], SCALEFACTOR10);
    x3 = L_shr_pos(im[s * 7], SCALEFACTOR10);
    x4 = L_shr_pos(im[s * 9], SCALEFACTOR10);

    s1  = L_add(x1, x4);
    s4  = L_sub(x1, x4);
    s3  = L_add(x3, x2);
    s2  = L_sub(x3, x2);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y03 = L_add(x0, s1);
    s1  = L_add(y03, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y06 = L_add(r1, s2);
    y18 = L_sub(r1, s2);
    y10 = L_sub(r3, s4);
    y14 = L_add(r3, s4);

    y07 = L_sub(s1, r2);
    y19 = L_add(s1, r2);
    y11 = L_add(s3, r4);
    y15 = L_sub(s3, r4);

    /* 5 fft2 stages */
    re[s * 0] = L_add(y00, y02); move32();
    im[s * 0] = L_add(y01, y03); move32();
    re[s * 5] = L_sub(y00, y02); move32();
    im[s * 5] = L_sub(y01, y03); move32();

    re[s * 2] = L_add(y04, y06); move32();
    im[s * 2] = L_add(y05, y07); move32();
    re[s * 7] = L_sub(y04, y06); move32();
    im[s * 7] = L_sub(y05, y07); move32();

    re[s * 4] = L_add(y08, y10); move32();
    im[s * 4] = L_add(y09, y11); move32();
    re[s * 9] = L_sub(y08, y10); move32();
    im[s * 9] = L_sub(y09, y11); move32();

    re[s * 6] = L_add(y12, y14); move32();
    im[s * 6] = L_add(y13, y15); move32();
    re[s * 1] = L_sub(y12, y14); move32();
    im[s * 1] = L_sub(y13, y15); move32();

    re[s * 8] = L_add(y16, y18); move32();
    im[s * 8] = L_add(y17, y19); move32();
    re[s * 3] = L_sub(y16, y18); move32();
    im[s * 3] = L_sub(y17, y19); move32();

    Dyn_Mem_Deluxe_Out();
}

/**
 * \brief    Function performs a complex 15-point FFT
 *           The FFT is performed inplace. The result of the FFT
 *           is scaled by SCALEFACTOR15 bits.
 *
 *           WOPS with 32x16 bit multiplications:  354 cycles
 *
 * \param    [i/o] re    real input / output
 * \param    [i/o] im    imag input / output
 * \param    [i  ] s     stride real and imag input / output
 *
 * \return   void
 */


static void fft15(Word32 *re, Word32 *im, Word16 s)
{
    Dyn_Mem_Deluxe_In(
        Word32 t;
        Word32 r1, r2, r3, r4;
        Word32 s1, s2, s3, s4;
        Word32 x00, x01, x02, x03, x04, x05, x06, x07, x08, x09;
        Word32 x10, x11, x12, x13, x14, x15, x16, x17, x18, x19;
        Word32 x20, x21, x22, x23, x24, x25, x26, x27, x28, x29;
        Word32 y00, y01, y02, y03, y04, y05, y06, y07, y08, y09;
        Word32 y10, y11, y12, y13, y14, y15, y16, y17, y18, y19;
        Word32 y20, y21, y22, y23, y24, y25, y26, y27, y28, y29;
    );

    x00 = L_shr_pos(re[s * 0], SCALEFACTOR15);
    x01 = L_shr_pos(im[s * 0], SCALEFACTOR15);
    x02 = L_shr_pos(re[s * 3], SCALEFACTOR15);
    x03 = L_shr_pos(im[s * 3], SCALEFACTOR15);
    x04 = L_shr_pos(re[s * 6], SCALEFACTOR15);
    x05 = L_shr_pos(im[s * 6], SCALEFACTOR15);
    x06 = L_shr_pos(re[s * 9], SCALEFACTOR15);
    x07 = L_shr_pos(im[s * 9], SCALEFACTOR15);
    x08 = L_shr_pos(re[s * 12], SCALEFACTOR15);
    x09 = L_shr_pos(im[s * 12], SCALEFACTOR15);

    x10 = L_shr_pos(re[s * 5], SCALEFACTOR15);
    x11 = L_shr_pos(im[s * 5], SCALEFACTOR15);
    x12 = L_shr_pos(re[s * 8], SCALEFACTOR15);
    x13 = L_shr_pos(im[s * 8], SCALEFACTOR15);
    x14 = L_shr_pos(re[s * 11], SCALEFACTOR15);
    x15 = L_shr_pos(im[s * 11], SCALEFACTOR15);
    x16 = L_shr_pos(re[s * 14], SCALEFACTOR15);
    x17 = L_shr_pos(im[s * 14], SCALEFACTOR15);
    x18 = L_shr_pos(re[s * 2], SCALEFACTOR15);
    x19 = L_shr_pos(im[s * 2], SCALEFACTOR15);

    x20 = L_shr_pos(re[s * 10], SCALEFACTOR15);
    x21 = L_shr_pos(im[s * 10], SCALEFACTOR15);
    x22 = L_shr_pos(re[s * 13], SCALEFACTOR15);
    x23 = L_shr_pos(im[s * 13], SCALEFACTOR15);
    x24 = L_shr_pos(re[s * 1], SCALEFACTOR15);
    x25 = L_shr_pos(im[s * 1], SCALEFACTOR15);
    x26 = L_shr_pos(re[s * 4], SCALEFACTOR15);
    x27 = L_shr_pos(im[s * 4], SCALEFACTOR15);
    x28 = L_shr_pos(re[s * 7], SCALEFACTOR15);
    x29 = L_shr_pos(im[s * 7], SCALEFACTOR15);

    /* 1. FFT5 stage */

    /* real part */
    r1  = L_add(x02, x08);
    r4  = L_sub(x02, x08);
    r3  = L_add(x04, x06);
    r2  = L_sub(x04, x06);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y00 = L_add(x00, r1);
    r1  = L_add(y00, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    s1  = L_add(x03, x09);
    s4  = L_sub(x03, x09);
    s3  = L_add(x05, x07);
    s2  = L_sub(x05, x07);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y01 = L_add(x01, s1);
    s1  = L_add(y01, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y02 = L_add(r1, s2);
    y08 = L_sub(r1, s2);
    y04 = L_sub(r3, s4);
    y06 = L_add(r3, s4);

    y03 = L_sub(s1, r2);
    y09 = L_add(s1, r2);
    y05 = L_add(s3, r4);
    y07 = L_sub(s3, r4);

    /* 2. FFT5 stage */

    /* real part */
    r1  = L_add(x12, x18);
    r4  = L_sub(x12, x18);
    r3  = L_add(x14, x16);
    r2  = L_sub(x14, x16);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y10 = L_add(x10, r1);
    r1  = L_add(y10, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    s1  = L_add(x13, x19);
    s4  = L_sub(x13, x19);
    s3  = L_add(x15, x17);
    s2  = L_sub(x15, x17);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y11 = L_add(x11, s1);
    s1  = L_add(y11, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y12 = L_add(r1, s2);
    y18 = L_sub(r1, s2);
    y14 = L_sub(r3, s4);
    y16 = L_add(r3, s4);

    y13 = L_sub(s1, r2);
    y19 = L_add(s1, r2);
    y15 = L_add(s3, r4);
    y17 = L_sub(s3, r4);

    /* 3. FFT5 stage */

    /* real part */
    r1  = L_add(x22, x28);
    r4  = L_sub(x22, x28);
    r3  = L_add(x24, x26);
    r2  = L_sub(x24, x26);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y20 = L_add(x20, r1);
    r1  = L_add(y20, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    s1  = L_add(x23, x29);
    s4  = L_sub(x23, x29);
    s3  = L_add(x25, x27);
    s2  = L_sub(x25, x27);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y21 = L_add(x21, s1);
    s1  = L_add(y21, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y22 = L_add(r1, s2);
    y28 = L_sub(r1, s2);
    y24 = L_sub(r3, s4);
    y26 = L_add(r3, s4);

    y23 = L_sub(s1, r2);
    y29 = L_add(s1, r2);
    y25 = L_add(s3, r4);
    y27 = L_sub(s3, r4);

    /* 1. FFT3 stage */

    /* real part */
    r1        = L_add(y10, y20);
    r2        = Mpy_32_xx(L_sub(y10, y20), C31);
    re[s * 0] = L_add(y00, r1); move32();
    r1        = L_sub(y00, L_shr_pos(r1, 1));

    /* imaginary part */
    s1        = L_add(y11, y21);
    s2        = Mpy_32_xx(L_sub(y11, y21), C31);
    im[s * 0] = L_add(y01, s1); move32();
    s1        = L_sub(y01, L_shr_pos(s1, 1));

    /* combination */
    re[s * 10] = L_sub(r1, s2); move32();
    re[s * 5]  = L_add(r1, s2); move32();
    im[s * 10] = L_add(s1, r2); move32();
    im[s * 5]  = L_sub(s1, r2); move32();

    /* 2. FFT3 stage */

    /* real part */
    r1        = L_add(y12, y22);
    r2        = Mpy_32_xx(L_sub(y12, y22), C31);
    re[s * 6] = L_add(y02, r1); move32();
    r1        = L_sub(y02, L_shr_pos(r1, 1));

    /* imaginary part */
    s1        = L_add(y13, y23);
    s2        = Mpy_32_xx(L_sub(y13, y23), C31);
    im[s * 6] = L_add(y03, s1); move32();
    s1        = L_sub(y03, L_shr_pos(s1, 1));

    /* combination */
    re[s * 1]  = L_sub(r1, s2); move32();
    re[s * 11] = L_add(r1, s2); move32();
    im[s * 1]  = L_add(s1, r2); move32();
    im[s * 11] = L_sub(s1, r2); move32();

    /* 3. FFT3 stage */

    /* real part */
    r1         = L_add(y14, y24);
    r2         = Mpy_32_xx(L_sub(y14, y24), C31);
    re[s * 12] = L_add(y04, r1); move32();
    r1         = L_sub(y04, L_shr_pos(r1, 1));

    /* imaginary part */
    s1         = L_add(y15, y25);
    s2         = Mpy_32_xx(L_sub(y15, y25), C31);
    im[s * 12] = L_add(y05, s1); move32();
    s1         = L_sub(y05, L_shr_pos(s1, 1));

    /* combination */
    re[s * 7] = L_sub(r1, s2); move32();
    re[s * 2] = L_add(r1, s2); move32();
    im[s * 7] = L_add(s1, r2); move32();
    im[s * 2] = L_sub(s1, r2); move32();

    /* 4. FFT3 stage */

    /* real part */
    r1        = L_add(y16, y26);
    r2        = Mpy_32_xx(L_sub(y16, y26), C31);
    re[s * 3] = L_add(y06, r1); move32();
    r1        = L_sub(y06, L_shr_pos(r1, 1));

    /* imaginary part */
    s1        = L_add(y17, y27);
    s2        = Mpy_32_xx(L_sub(y17, y27), C31);
    im[s * 3] = L_add(y07, s1); move32();
    s1        = L_sub(y07, L_shr_pos(s1, 1));

    /* combination */
    re[s * 13] = L_sub(r1, s2); move32();
    re[s * 8]  = L_add(r1, s2); move32();
    im[s * 13] = L_add(s1, r2); move32();
    im[s * 8]  = L_sub(s1, r2); move32();

    /* 5. FFT3 stage */

    /* real part */
    r1        = L_add(y18, y28);
    r2        = Mpy_32_xx(L_sub(y18, y28), C31);
    re[s * 9] = L_add(y08, r1); move32();
    r1        = L_sub(y08, L_shr_pos(r1, 1));

    /* imaginary part */
    s1        = L_add(y19, y29);
    s2        = Mpy_32_xx(L_sub(y19, y29), C31);
    im[s * 9] = L_add(y09, s1); move32();
    s1        = L_sub(y09, L_shr_pos(s1, 1));

    /* combination */
    re[s * 4]  = L_sub(r1, s2); move32();
    re[s * 14] = L_add(r1, s2); move32();
    im[s * 4]  = L_add(s1, r2); move32();
    im[s * 14] = L_sub(s1, r2); move32();

    Dyn_Mem_Deluxe_Out();
}


static void fft12(Word32 *pInput)
{
    Dyn_Mem_Deluxe_In(
        Word32  aDst[24];
        Word32 *pSrc, *pDst;
        Counter i;
        Word32  r1, r2, s1, s2, pD;
        Word32  re, im;
        Word16  vre, vim;
    );

    pSrc = pInput; move16();
    pDst = aDst;   move16();

    /* First 3*2 samples are shifted right by 2 before output */
    r1      = L_add(L_shr_pos(pSrc[8], 2), L_shr_pos(pSrc[16], 2));
    r2      = Mpy_32_16(L_sub(L_shr_pos(pSrc[8], 2), L_shr_pos(pSrc[16], 2)), C31);
    pD      = L_shr_pos(pSrc[0], 2);
    pDst[0] = L_shr_pos(L_add(pD, r1), 1);
    r1      = L_sub(pD, L_shr_pos(r1, 1));

    /* imaginary part */
    s1      = L_add(L_shr_pos(pSrc[9], 2), L_shr_pos(pSrc[17], 2));
    s2      = Mpy_32_16(L_sub(L_shr_pos(pSrc[9], 2), L_shr_pos(pSrc[17], 2)), C31);
    pD      = L_shr_pos(pSrc[1], 2);
    pDst[1] = L_shr_pos(L_add(pD, s1), 1);
    s1      = L_sub(pD, L_shr_pos(s1, 1));

    r1 = L_shr_pos(r1, 1);
    r2 = L_shr_pos(r2, 1);
    s1 = L_shr_pos(s1, 1);
    s2 = L_shr_pos(s2, 1);

    /* combination */
    pDst[2] = L_sub(r1, s2);
    pDst[3] = L_add(s1, r2);
    pDst[4] = L_add(r1, s2);
    pDst[5] = L_sub(s1, r2);
    pSrc += 2;
    pDst += 6;

    vre = add(0x6eda, 0);
    vim = add(0x4000, 0);

    FOR (i = 0; i < 2; i++)
    {
        /* sample 0,1 are shifted right by 2 before output */
        /* sample 2,3 4,5 are shifted right by 1 and complex multiplied before output */

        r1      = L_add(L_shr_pos(pSrc[8], 2), L_shr_pos(pSrc[16], 2));
        r2      = Mpy_32_16(L_sub(L_shr_pos(pSrc[8], 2), L_shr_pos(pSrc[16], 2)), C31);
        pD      = L_shr_pos(pSrc[0], 2);
        pDst[0] = L_shr_pos(L_add(pD, r1), 1);
        r1      = L_sub(pD, L_shr_pos(r1, 1));

        /* imaginary part */
        s1      = L_add(L_shr_pos(pSrc[9], 2), L_shr_pos(pSrc[17], 2));
        s2      = Mpy_32_16(L_sub(L_shr_pos(pSrc[9], 2), L_shr_pos(pSrc[17], 2)), C31);
        pD      = L_shr_pos(pSrc[1], 2);
        pDst[1] = L_shr_pos(L_add(pD, s1), 1);
        s1      = L_sub(pD, L_shr_pos(s1, 1));

        r1 = L_shr_pos(r1, 1);
        r2 = L_shr_pos(r2, 1);
        s1 = L_shr_pos(s1, 1);
        s2 = L_shr_pos(s2, 1);

        /* combination */
        re = L_sub(r1, s2);
        im = L_add(s1, r2);
        cplxMpy_32_16(&pDst[3], &pDst[2], im, re, vre, vim);
        re  = L_add(r1, s2);
        im  = L_sub(s1, r2);
        vre = add(0x4000, 0);
        if (i == 1)
            vre = negate(vre); /* 0xC000 */
        if (i == 0)
            vim = add(0x6eda, 0);
        cplxMpy_32_16(&pDst[5], &pDst[4], im, re, vre, vim);

        pDst += 6;
        pSrc += 2;
    }
    /* sample 0,1 are shifted right by 2 before output */
    /* sample 2,3 is shifted right by 1 and complex multiplied with (0.0,+1.0) */
    /* sample 4,5 is shifted right by 1 and complex multiplied with (-1.0,0.0) */
    r1      = L_add(L_shr_pos(pSrc[8], 2), L_shr_pos(pSrc[16], 2));
    r2      = Mpy_32_16(L_sub(L_shr_pos(pSrc[8], 2), L_shr_pos(pSrc[16], 2)), C31);
    pD      = L_shr_pos(pSrc[0], 2);
    pDst[0] = L_shr_pos(L_add(pD, r1), 1);
    r1      = L_sub(pD, L_shr_pos(r1, 1));

    /* imaginary part */
    s1      = L_add(L_shr_pos(pSrc[9], 2), L_shr_pos(pSrc[17], 2));
    s2      = Mpy_32_16(L_sub(L_shr_pos(pSrc[9], 2), L_shr_pos(pSrc[17], 2)), C31);
    pD      = L_shr_pos(pSrc[1], 2);
    pDst[1] = L_shr_pos(L_add(pD, s1), 1);
    s1      = L_sub(pD, L_shr_pos(s1, 1));

    r1 = L_shr_pos(r1, 1);
    r2 = L_shr_pos(r2, 1);
    s1 = L_shr_pos(s1, 1);
    s2 = L_shr_pos(s2, 1);

    /* combination */

    pDst[2] = L_add(s1, r2);           move32();
    pDst[3] = L_sub(s2, r1);           move32();
    pDst[4] = L_negate(L_add(r1, s2)); move32();
    pDst[5] = L_sub(r2, s1);           move32();
    /* Perform 3 times the fft of length 4. The input samples are at the address of aDst and the
     output samples are at the address of pInput. The input vector for the fft of length 4 is built
     of the interleaved samples in aDst, the output samples are stored consecutively at the address
     of pInput.
     */
    move16(); move16();
    pSrc = aDst;
    pDst = pInput;
    FOR (i = 0; i < 3; i++)
    {
        /* inline FFT4 merged with incoming resorting loop */
        r1 = L_add(L_shr_pos(pSrc[0], 2), L_shr_pos(pSrc[12], 2)); /* Re A + Re B */
        r2 = L_add(L_shr_pos(pSrc[6], 2), L_shr_pos(pSrc[18], 2)); /* Re C + Re D */
        s1 = L_add(L_shr_pos(pSrc[1], 2), L_shr_pos(pSrc[13], 2)); /* Im A + Im B */
        s2 = L_add(L_shr_pos(pSrc[7], 2), L_shr_pos(pSrc[19], 2)); /* Im C + Im D */

        pDst[0] = L_add(r1, r2); /* Re A' = Re A + Re B + Re C + Re D */
        pDst[1] = L_add(s1, s2); /* Im A' = Im A + Im B + Im C + Im D */

        re = L_sub(r1, L_shr_pos(pSrc[12], 1)); /* Re A - Re B */
        im = L_sub(s1, L_shr_pos(pSrc[13], 1)); /* Im A - Im B */

        pDst[12] = L_sub(r1, r2); /* Re C' = Re A + Re B - Re C - Re D */
        pDst[13] = L_sub(s1, s2); /* Im C' = Im A + Im B - Im C - Im D */

        r2 = L_sub(r2, L_shr_pos(pSrc[18], 1)); /* Re C - Re D */
        s2 = L_sub(s2, L_shr_pos(pSrc[19], 1)); /* Im C - Im D */

        pDst[6]  = L_add(re, s2); /* Re B' = Re A - Re B + Im C - Im D */
        pDst[18] = L_sub(re, s2); /* Re D' = Re A - Re B - Im C + Im D */
        pDst[7]  = L_sub(im, r2); /* Im B' = Im A - Im B - Re C + Re D */
        pDst[19] = L_add(im, r2); /* Im D' = Im A - Im B + Re C - Re D */

        pSrc += 2;
        pDst += 2;
    }

    Dyn_Mem_Deluxe_Out();
}

/**
 * \brief    Function performs a complex 16-point FFT
 *           The FFT is performed inplace. The result of the FFT
 *           is scaled by SCALEFACTOR16 bits.
 *
 *           WOPS with 32x16 bit multiplications (scale on ):  288 cycles
 *           WOPS with 32x16 bit multiplications (scale off):  256 cycles
 *
 * \param    [i/o] re    real input / output
 * \param    [i/o] im    imag input / output
 * \param    [i  ] s     stride real and imag input / output
 *
 * \return   void
 */


static void fft16(Word32 *re, Word32 *im, Word16 s)
{
    Dyn_Mem_Deluxe_In(
        Word32 x0, x1, x2, x3, x4, x5, x6, x7;
        Word32 t0, t1, t2, t3, t4, t5, t6, t7;
        Word32 y00, y01, y02, y03, y04, y05, y06, y07;
        Word32 y08, y09, y10, y11, y12, y13, y14, y15;
        Word32 y16, y17, y18, y19, y20, y21, y22, y23;
        Word32 y24, y25, y26, y27, y28, y29, y30, y31;
    );

    x0 = L_shr_pos(re[s * 0], SCALEFACTOR16);
    x1 = L_shr_pos(im[s * 0], SCALEFACTOR16);
    x2 = L_shr_pos(re[s * 4], SCALEFACTOR16);
    x3 = L_shr_pos(im[s * 4], SCALEFACTOR16);
    x4 = L_shr_pos(re[s * 8], SCALEFACTOR16);
    x5 = L_shr_pos(im[s * 8], SCALEFACTOR16);
    x6 = L_shr_pos(re[s * 12], SCALEFACTOR16);
    x7 = L_shr_pos(im[s * 12], SCALEFACTOR16);

    /* Pre-additions */
    t0 = L_add(x0, x4);
    t2 = L_sub(x0, x4);
    t1 = L_add(x1, x5);
    t3 = L_sub(x1, x5);
    t4 = L_add(x2, x6);
    t7 = L_sub(x2, x6);
    t5 = L_add(x7, x3);
    t6 = L_sub(x7, x3);

    /* Post-additions */
    y00 = L_add(t0, t4);
    y01 = L_add(t1, t5);
    y02 = L_sub(t2, t6);
    y03 = L_sub(t3, t7);
    y04 = L_sub(t0, t4);
    y05 = L_sub(t1, t5);
    y06 = L_add(t2, t6);
    y07 = L_add(t3, t7);

    x0 = L_shr_pos(re[s * 1], SCALEFACTOR16);
    x1 = L_shr_pos(im[s * 1], SCALEFACTOR16);
    x2 = L_shr_pos(re[s * 5], SCALEFACTOR16);
    x3 = L_shr_pos(im[s * 5], SCALEFACTOR16);
    x4 = L_shr_pos(re[s * 9], SCALEFACTOR16);
    x5 = L_shr_pos(im[s * 9], SCALEFACTOR16);
    x6 = L_shr_pos(re[s * 13], SCALEFACTOR16);
    x7 = L_shr_pos(im[s * 13], SCALEFACTOR16);

    /* Pre-additions */
    t0 = L_add(x0, x4);
    t2 = L_sub(x0, x4);
    t1 = L_add(x1, x5);
    t3 = L_sub(x1, x5);
    t4 = L_add(x2, x6);
    t7 = L_sub(x2, x6);
    t5 = L_add(x7, x3);
    t6 = L_sub(x7, x3);

    /* Post-additions */
    y08 = L_add(t0, t4);
    y09 = L_add(t1, t5);
    y10 = L_sub(t2, t6);
    y11 = L_sub(t3, t7);
    y12 = L_sub(t0, t4);
    y13 = L_sub(t1, t5);
    y14 = L_add(t2, t6);
    y15 = L_add(t3, t7);

    x0 = L_shr_pos(re[s * 2], SCALEFACTOR16);
    x1 = L_shr_pos(im[s * 2], SCALEFACTOR16);
    x2 = L_shr_pos(re[s * 6], SCALEFACTOR16);
    x3 = L_shr_pos(im[s * 6], SCALEFACTOR16);
    x4 = L_shr_pos(re[s * 10], SCALEFACTOR16);
    x5 = L_shr_pos(im[s * 10], SCALEFACTOR16);
    x6 = L_shr_pos(re[s * 14], SCALEFACTOR16);
    x7 = L_shr_pos(im[s * 14], SCALEFACTOR16);

    /* Pre-additions */
    t0 = L_add(x0, x4);
    t2 = L_sub(x0, x4);
    t1 = L_add(x1, x5);
    t3 = L_sub(x1, x5);
    t4 = L_add(x2, x6);
    t7 = L_sub(x2, x6);
    t5 = L_add(x7, x3);
    t6 = L_sub(x7, x3);

    /* Post-additions */
    y16 = L_add(t0, t4);
    y17 = L_add(t1, t5);
    y18 = L_sub(t2, t6);
    y19 = L_sub(t3, t7);
    y20 = L_sub(t1, t5);
    y21 = L_sub(t4, t0);
    y22 = L_add(t2, t6);
    y23 = L_add(t3, t7);

    x0 = L_shr_pos(re[s * 3], SCALEFACTOR16);
    x1 = L_shr_pos(im[s * 3], SCALEFACTOR16);
    x2 = L_shr_pos(re[s * 7], SCALEFACTOR16);
    x3 = L_shr_pos(im[s * 7], SCALEFACTOR16);
    x4 = L_shr_pos(re[s * 11], SCALEFACTOR16);
    x5 = L_shr_pos(im[s * 11], SCALEFACTOR16);
    x6 = L_shr_pos(re[s * 15], SCALEFACTOR16);
    x7 = L_shr_pos(im[s * 15], SCALEFACTOR16);

    /* Pre-additions */
    t0 = L_add(x0, x4);
    t2 = L_sub(x0, x4);
    t1 = L_add(x1, x5);
    t3 = L_sub(x1, x5);
    t4 = L_add(x2, x6);
    t7 = L_sub(x2, x6);
    t5 = L_add(x7, x3);
    t6 = L_sub(x7, x3);

    /* Post-additions */
    y24 = L_add(t0, t4);
    y25 = L_add(t1, t5);
    y26 = L_sub(t2, t6);
    y27 = L_sub(t3, t7);
    y28 = L_sub(t0, t4);
    y29 = L_sub(t1, t5);
    y30 = L_add(t2, t6);
    y31 = L_add(t3, t7);

    /* rotation */

    x0  = Mpy_32_xx(y22, C162);
    x1  = Mpy_32_xx(y23, C162);
    y22 = L_sub(x0, x1);
    y23 = L_add(x0, x1);

    x0  = Mpy_32_xx(y28, C162);
    x1  = Mpy_32_xx(y29, C162);
    y28 = L_sub(x0, x1);
    y29 = L_add(x0, x1);

    x0  = Mpy_32_xx(y12, C161);
    x1  = Mpy_32_xx(y13, C161);
    y12 = L_add(x0, x1);
    y13 = L_sub(x1, x0);

    x0  = Mpy_32_xx(y18, C161);
    x1  = Mpy_32_xx(y19, C161);
    y18 = L_add(x0, x1);
    y19 = L_sub(x1, x0);

    x0  = Mpy_32_xx(y10, C163);
    x1  = Mpy_32_xx(y11, C166);
    x2  = Mpy_32_xx(y10, C166);
    x3  = Mpy_32_xx(y11, C163);
    y10 = L_sub(x0, x1);
    y11 = L_add(x2, x3);

    x0  = Mpy_32_xx(y14, C165);
    x1  = Mpy_32_xx(y15, C164);
    x2  = Mpy_32_xx(y14, C164);
    x3  = Mpy_32_xx(y15, C165);
    y14 = L_sub(x0, x1);
    y15 = L_add(x2, x3);

    x0  = Mpy_32_xx(y26, C165);
    x1  = Mpy_32_xx(y27, C164);
    x2  = Mpy_32_xx(y26, C164);
    x3  = Mpy_32_xx(y27, C165);
    y26 = L_sub(x0, x1);
    y27 = L_add(x2, x3);

    x0  = Mpy_32_xx(y30, C164);
    x1  = Mpy_32_xx(y31, C165);
    x2  = Mpy_32_xx(y30, C165);
    x3  = Mpy_32_xx(y31, C164);
    y30 = L_sub(x0, x1);
    y31 = L_add(x2, x3);

    /* Pre-additions */

    t0 = L_add(y00, y16);
    t2 = L_sub(y00, y16);
    t1 = L_add(y01, y17);
    t3 = L_sub(y01, y17);
    t4 = L_add(y08, y24);
    t7 = L_sub(y08, y24);
    t5 = L_add(y25, y09);
    t6 = L_sub(y25, y09);

    /* Post-additions */

    re[s * 0]  = L_add(t0, t4); move32();
    im[s * 0]  = L_add(t1, t5); move32();
    re[s * 4]  = L_sub(t2, t6); move32();
    im[s * 4]  = L_sub(t3, t7); move32();
    re[s * 8]  = L_sub(t0, t4); move32();
    im[s * 8]  = L_sub(t1, t5); move32();
    re[s * 12] = L_add(t2, t6); move32();
    im[s * 12] = L_add(t3, t7); move32();

    /* Pre-additions */

    t0 = L_add(y02, y18);
    t2 = L_sub(y02, y18);
    t1 = L_add(y03, y19);
    t3 = L_sub(y03, y19);
    t4 = L_add(y10, y26);
    t7 = L_sub(y10, y26);
    t5 = L_add(y27, y11);
    t6 = L_sub(y27, y11);

    /* Post-additions */

    re[s * 1]  = L_add(t0, t4); move32();
    im[s * 1]  = L_add(t1, t5); move32();
    re[s * 5]  = L_sub(t2, t6); move32();
    im[s * 5]  = L_sub(t3, t7); move32();
    re[s * 9]  = L_sub(t0, t4); move32();
    im[s * 9]  = L_sub(t1, t5); move32();
    re[s * 13] = L_add(t2, t6); move32();
    im[s * 13] = L_add(t3, t7); move32();

    /* Pre-additions */

    t0 = L_add(y04, y20);
    t2 = L_sub(y04, y20);
    t1 = L_add(y05, y21);
    t3 = L_sub(y05, y21);
    t4 = L_add(y12, y28);
    t7 = L_sub(y12, y28);
    t5 = L_add(y29, y13);
    t6 = L_sub(y29, y13);

    /* Post-additions */

    re[s * 2]  = L_add(t0, t4); move32();
    im[s * 2]  = L_add(t1, t5); move32();
    re[s * 6]  = L_sub(t2, t6); move32();
    im[s * 6]  = L_sub(t3, t7); move32();
    re[s * 10] = L_sub(t0, t4); move32();
    im[s * 10] = L_sub(t1, t5); move32();
    re[s * 14] = L_add(t2, t6); move32();
    im[s * 14] = L_add(t3, t7); move32();

    /* Pre-additions */

    t0 = L_add(y06, y22);
    t2 = L_sub(y06, y22);
    t1 = L_add(y07, y23);
    t3 = L_sub(y07, y23);
    t4 = L_add(y14, y30);
    t7 = L_sub(y14, y30);
    t5 = L_add(y31, y15);
    t6 = L_sub(y31, y15);

    /* Post-additions */

    re[s * 3]  = L_add(t0, t4); move32();
    im[s * 3]  = L_add(t1, t5); move32();
    re[s * 7]  = L_sub(t2, t6); move32();
    im[s * 7]  = L_sub(t3, t7); move32();
    re[s * 11] = L_sub(t0, t4); move32();
    im[s * 11] = L_sub(t1, t5); move32();
    re[s * 15] = L_add(t2, t6); move32();
    im[s * 15] = L_add(t3, t7); move32();

    Dyn_Mem_Deluxe_Out();
}

/**
 * \brief    Function performs a complex 20-point FFT
 *           The FFT is performed inplace. The result of the FFT
 *           is scaled by SCALEFACTOR20 bits.
 *
 *           WOPS with 32x16 bit multiplications:  432 cycles
 *
 * \param    [i/o] re    real input / output
 * \param    [i/o] im    imag input / output
 * \param    [i  ] s     stride real and imag input / output
 *
 * \return   void
 */


static void fft20(Word32 *re, Word32 *im, Word16 s)
{
    Dyn_Mem_Deluxe_In(
        Word32 r1, r2, r3, r4;
        Word32 s1, s2, s3, s4;
        Word32 x0, x1, x2, x3, x4;
        Word32 t, t0, t1, t2, t3, t4, t5, t6, t7;
        Word32 y00, y01, y02, y03, y04, y05, y06, y07, y08, y09;
        Word32 y10, y11, y12, y13, y14, y15, y16, y17, y18, y19;
        Word32 y20, y21, y22, y23, y24, y25, y26, y27, y28, y29;
        Word32 y30, y31, y32, y33, y34, y35, y36, y37, y38, y39;
    );

    /* 1. FFT5 stage */

    /* real part */
    x0 = L_shr_pos(re[s * 0], SCALEFACTOR20);
    x1 = L_shr_pos(re[s * 16], SCALEFACTOR20);
    x2 = L_shr_pos(re[s * 12], SCALEFACTOR20);
    x3 = L_shr_pos(re[s * 8], SCALEFACTOR20);
    x4 = L_shr_pos(re[s * 4], SCALEFACTOR20);

    r1  = L_add(x1, x4);
    r4  = L_sub(x1, x4);
    r3  = L_add(x2, x3);
    r2  = L_sub(x2, x3);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y00 = L_add(x0, r1);
    r1  = L_add(y00, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    x0 = L_shr_pos(im[s * 0], SCALEFACTOR20);
    x1 = L_shr_pos(im[s * 16], SCALEFACTOR20);
    x2 = L_shr_pos(im[s * 12], SCALEFACTOR20);
    x3 = L_shr_pos(im[s * 8], SCALEFACTOR20);
    x4 = L_shr_pos(im[s * 4], SCALEFACTOR20);

    s1  = L_add(x1, x4);
    s4  = L_sub(x1, x4);
    s3  = L_add(x2, x3);
    s2  = L_sub(x2, x3);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y01 = L_add(x0, s1);
    s1  = L_add(y01, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y08 = L_add(r1, s2);
    y32 = L_sub(r1, s2);
    y16 = L_sub(r3, s4);
    y24 = L_add(r3, s4);

    y09 = L_sub(s1, r2);
    y33 = L_add(s1, r2);
    y17 = L_add(s3, r4);
    y25 = L_sub(s3, r4);

    /* 2. FFT5 stage */

    /* real part */
    x0 = L_shr_pos(re[s * 5], SCALEFACTOR20);
    x1 = L_shr_pos(re[s * 1], SCALEFACTOR20);
    x2 = L_shr_pos(re[s * 17], SCALEFACTOR20);
    x3 = L_shr_pos(re[s * 13], SCALEFACTOR20);
    x4 = L_shr_pos(re[s * 9], SCALEFACTOR20);

    r1  = L_add(x1, x4);
    r4  = L_sub(x1, x4);
    r3  = L_add(x2, x3);
    r2  = L_sub(x2, x3);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y02 = L_add(x0, r1);
    r1  = L_add(y02, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    x0 = L_shr_pos(im[s * 5], SCALEFACTOR20);
    x1 = L_shr_pos(im[s * 1], SCALEFACTOR20);
    x2 = L_shr_pos(im[s * 17], SCALEFACTOR20);
    x3 = L_shr_pos(im[s * 13], SCALEFACTOR20);
    x4 = L_shr_pos(im[s * 9], SCALEFACTOR20);

    s1  = L_add(x1, x4);
    s4  = L_sub(x1, x4);
    s3  = L_add(x2, x3);
    s2  = L_sub(x2, x3);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y03 = L_add(x0, s1);
    s1  = L_add(y03, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y10 = L_add(r1, s2);
    y34 = L_sub(r1, s2);
    y18 = L_sub(r3, s4);
    y26 = L_add(r3, s4);

    y11 = L_sub(s1, r2);
    y35 = L_add(s1, r2);
    y19 = L_add(s3, r4);
    y27 = L_sub(s3, r4);

    /* 3. FFT5 stage */

    /* real part */
    x0 = L_shr_pos(re[s * 10], SCALEFACTOR20);
    x1 = L_shr_pos(re[s * 6], SCALEFACTOR20);
    x2 = L_shr_pos(re[s * 2], SCALEFACTOR20);
    x3 = L_shr_pos(re[s * 18], SCALEFACTOR20);
    x4 = L_shr_pos(re[s * 14], SCALEFACTOR20);

    r1  = L_add(x1, x4);
    r4  = L_sub(x1, x4);
    r3  = L_add(x2, x3);
    r2  = L_sub(x2, x3);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y04 = L_add(x0, r1);
    r1  = L_add(y04, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    x0 = L_shr_pos(im[s * 10], SCALEFACTOR20);
    x1 = L_shr_pos(im[s * 6], SCALEFACTOR20);
    x2 = L_shr_pos(im[s * 2], SCALEFACTOR20);
    x3 = L_shr_pos(im[s * 18], SCALEFACTOR20);
    x4 = L_shr_pos(im[s * 14], SCALEFACTOR20);

    s1  = L_add(x1, x4);
    s4  = L_sub(x1, x4);
    s3  = L_add(x2, x3);
    s2  = L_sub(x2, x3);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y05 = L_add(x0, s1);
    s1  = L_add(y05, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y12 = L_add(r1, s2);
    y36 = L_sub(r1, s2);
    y20 = L_sub(r3, s4);
    y28 = L_add(r3, s4);

    y13 = L_sub(s1, r2);
    y37 = L_add(s1, r2);
    y21 = L_add(s3, r4);
    y29 = L_sub(s3, r4);

    /* 4. FFT5 stage */

    /* real part */
    x0 = L_shr_pos(re[s * 15], SCALEFACTOR20);
    x1 = L_shr_pos(re[s * 11], SCALEFACTOR20);
    x2 = L_shr_pos(re[s * 7], SCALEFACTOR20);
    x3 = L_shr_pos(re[s * 3], SCALEFACTOR20);
    x4 = L_shr_pos(re[s * 19], SCALEFACTOR20);

    r1  = L_add(x1, x4);
    r4  = L_sub(x1, x4);
    r3  = L_add(x2, x3);
    r2  = L_sub(x2, x3);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y06 = L_add(x0, r1);
    r1  = L_add(y06, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    x0 = L_shr_pos(im[s * 15], SCALEFACTOR20);
    x1 = L_shr_pos(im[s * 11], SCALEFACTOR20);
    x2 = L_shr_pos(im[s * 7], SCALEFACTOR20);
    x3 = L_shr_pos(im[s * 3], SCALEFACTOR20);
    x4 = L_shr_pos(im[s * 19], SCALEFACTOR20);

    s1  = L_add(x1, x4);
    s4  = L_sub(x1, x4);
    s3  = L_add(x2, x3);
    s2  = L_sub(x2, x3);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y07 = L_add(x0, s1);
    s1  = L_add(y07, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y14 = L_add(r1, s2);
    y38 = L_sub(r1, s2);
    y22 = L_sub(r3, s4);
    y30 = L_add(r3, s4);

    y15 = L_sub(s1, r2);
    y39 = L_add(s1, r2);
    y23 = L_add(s3, r4);
    y31 = L_sub(s3, r4);

    /* 1. FFT4 stage */

    /* Pre-additions */
    t0 = L_add(y00, y04);
    t2 = L_sub(y00, y04);
    t1 = L_add(y01, y05);
    t3 = L_sub(y01, y05);
    t4 = L_add(y02, y06);
    t7 = L_sub(y02, y06);
    t5 = L_add(y07, y03);
    t6 = L_sub(y07, y03);

    /* Post-additions */
    re[s * 0]  = L_add(t0, t4); move32();
    im[s * 0]  = L_add(t1, t5); move32();
    re[s * 5]  = L_sub(t2, t6); move32();
    im[s * 5]  = L_sub(t3, t7); move32();
    re[s * 10] = L_sub(t0, t4); move32();
    im[s * 10] = L_sub(t1, t5); move32();
    re[s * 15] = L_add(t2, t6); move32();
    im[s * 15] = L_add(t3, t7); move32();

    /* 2. FFT4 stage */

    /* Pre-additions */
    t0 = L_add(y08, y12);
    t2 = L_sub(y08, y12);
    t1 = L_add(y09, y13);
    t3 = L_sub(y09, y13);
    t4 = L_add(y10, y14);
    t7 = L_sub(y10, y14);
    t5 = L_add(y15, y11);
    t6 = L_sub(y15, y11);

    /* Post-additions */
    re[s * 4]  = L_add(t0, t4); move32();
    im[s * 4]  = L_add(t1, t5); move32();
    re[s * 9]  = L_sub(t2, t6); move32();
    im[s * 9]  = L_sub(t3, t7); move32();
    re[s * 14] = L_sub(t0, t4); move32();
    im[s * 14] = L_sub(t1, t5); move32();
    re[s * 19] = L_add(t2, t6); move32();
    im[s * 19] = L_add(t3, t7); move32();

    /* 3. FFT4 stage */

    /* Pre-additions */
    t0 = L_add(y16, y20);
    t2 = L_sub(y16, y20);
    t1 = L_add(y17, y21);
    t3 = L_sub(y17, y21);
    t4 = L_add(y18, y22);
    t7 = L_sub(y18, y22);
    t5 = L_add(y23, y19);
    t6 = L_sub(y23, y19);

    /* Post-additions */
    re[s * 8]  = L_add(t0, t4); move32();
    im[s * 8]  = L_add(t1, t5); move32();
    re[s * 13] = L_sub(t2, t6); move32();
    im[s * 13] = L_sub(t3, t7); move32();
    re[s * 18] = L_sub(t0, t4); move32();
    im[s * 18] = L_sub(t1, t5); move32();
    re[s * 3]  = L_add(t2, t6); move32();
    im[s * 3]  = L_add(t3, t7); move32();

    /* 4. FFT4 stage */

    /* Pre-additions */
    t0 = L_add(y24, y28);
    t2 = L_sub(y24, y28);
    t1 = L_add(y25, y29);
    t3 = L_sub(y25, y29);
    t4 = L_add(y26, y30);
    t7 = L_sub(y26, y30);
    t5 = L_add(y31, y27);
    t6 = L_sub(y31, y27);

    /* Post-additions */
    re[s * 12] = L_add(t0, t4); move32();
    im[s * 12] = L_add(t1, t5); move32();
    re[s * 17] = L_sub(t2, t6); move32();
    im[s * 17] = L_sub(t3, t7); move32();
    re[s * 2]  = L_sub(t0, t4); move32();
    im[s * 2]  = L_sub(t1, t5); move32();
    re[s * 7]  = L_add(t2, t6); move32();
    im[s * 7]  = L_add(t3, t7); move32();

    /* 5. FFT4 stage */

    /* Pre-additions */
    t0 = L_add(y32, y36);
    t2 = L_sub(y32, y36);
    t1 = L_add(y33, y37);
    t3 = L_sub(y33, y37);
    t4 = L_add(y34, y38);
    t7 = L_sub(y34, y38);
    t5 = L_add(y39, y35);
    t6 = L_sub(y39, y35);

    /* Post-additions */
    re[s * 16] = L_add(t0, t4); move32();
    im[s * 16] = L_add(t1, t5); move32();
    re[s * 1]  = L_sub(t2, t6); move32();
    im[s * 1]  = L_sub(t3, t7); move32();
    re[s * 6]  = L_sub(t0, t4); move32();
    im[s * 6]  = L_sub(t1, t5); move32();
    re[s * 11] = L_add(t2, t6); move32();
    im[s * 11] = L_add(t3, t7); move32();

    Dyn_Mem_Deluxe_Out();
}

/**
 * \brief    Function performs a complex 30-point FFT
 *           The FFT is performed inplace. The result of the FFT
 *           is scaled by SCALEFACTOR30 bits.
 *
 *           WOPS with 32x16 bit multiplications:  828 cycles
 *
 * \param    [i/o] re    real input / output
 * \param    [i/o] im    imag input / output
 * \param    [i  ] s     stride real and imag input / output
 *
 * \return   void
 */


static void fft30(Word32 *re, Word32 *im, Word16 s)
{
    Dyn_Mem_Deluxe_In(
        Word32 t;
        Word32 r1, r2, r3, r4;
        Word32 s1, s2, s3, s4;
        Word32 x00, x01, x02, x03, x04, x05, x06, x07, x08, x09;
        Word32 x10, x11, x12, x13, x14, x15, x16, x17, x18, x19;
        Word32 x20, x21, x22, x23, x24, x25, x26, x27, x28, x29;

        Word32 y00, y01, y02, y03, y04, y05, y06, y07, y08, y09;
        Word32 y10, y11, y12, y13, y14, y15, y16, y17, y18, y19;
        Word32 y20, y21, y22, y23, y24, y25, y26, y27, y28, y29;

        Word32 z00, z01, z02, z03, z04, z05, z06, z07, z08, z09;
        Word32 z10, z11, z12, z13, z14, z15, z16, z17, z18, z19;
        Word32 z20, z21, z22, z23, z24, z25, z26, z27, z28, z29;
        Word32 z30, z31, z32, z33, z34, z35, z36, z37, z38, z39;
        Word32 z40, z41, z42, z43, z44, z45, z46, z47, z48, z49;
        Word32 z50, z51, z52, z53, z54, z55, z56, z57, z58, z59;

        Word32 *rel, *reh, *iml, *imh;
    );

    rel = &re[s * 0];
    reh = &re[s * 15];
    iml = &im[s * 0];
    imh = &im[s * 15];

    /* 1. FFT15 stage */
    x00 = L_shr_pos(re[s * 0], SCALEFACTOR30_1);
    x01 = L_shr_pos(im[s * 0], SCALEFACTOR30_1);
    x02 = L_shr_pos(re[s * 18], SCALEFACTOR30_1);
    x03 = L_shr_pos(im[s * 18], SCALEFACTOR30_1);
    x04 = L_shr_pos(re[s * 6], SCALEFACTOR30_1);
    x05 = L_shr_pos(im[s * 6], SCALEFACTOR30_1);
    x06 = L_shr_pos(re[s * 24], SCALEFACTOR30_1);
    x07 = L_shr_pos(im[s * 24], SCALEFACTOR30_1);
    x08 = L_shr_pos(re[s * 12], SCALEFACTOR30_1);
    x09 = L_shr_pos(im[s * 12], SCALEFACTOR30_1);

    x10 = L_shr_pos(re[s * 20], SCALEFACTOR30_1);
    x11 = L_shr_pos(im[s * 20], SCALEFACTOR30_1);
    x12 = L_shr_pos(re[s * 8], SCALEFACTOR30_1);
    x13 = L_shr_pos(im[s * 8], SCALEFACTOR30_1);
    x14 = L_shr_pos(re[s * 26], SCALEFACTOR30_1);
    x15 = L_shr_pos(im[s * 26], SCALEFACTOR30_1);
    x16 = L_shr_pos(re[s * 14], SCALEFACTOR30_1);
    x17 = L_shr_pos(im[s * 14], SCALEFACTOR30_1);
    x18 = L_shr_pos(re[s * 2], SCALEFACTOR30_1);
    x19 = L_shr_pos(im[s * 2], SCALEFACTOR30_1);

    x20 = L_shr_pos(re[s * 10], SCALEFACTOR30_1);
    x21 = L_shr_pos(im[s * 10], SCALEFACTOR30_1);
    x22 = L_shr_pos(re[s * 28], SCALEFACTOR30_1);
    x23 = L_shr_pos(im[s * 28], SCALEFACTOR30_1);
    x24 = L_shr_pos(re[s * 16], SCALEFACTOR30_1);
    x25 = L_shr_pos(im[s * 16], SCALEFACTOR30_1);
    x26 = L_shr_pos(re[s * 4], SCALEFACTOR30_1);
    x27 = L_shr_pos(im[s * 4], SCALEFACTOR30_1);
    x28 = L_shr_pos(re[s * 22], SCALEFACTOR30_1);
    x29 = L_shr_pos(im[s * 22], SCALEFACTOR30_1);

    /* 1. FFT5 stage */

    /* real part */
    r1  = L_add(x02, x08);
    r4  = L_sub(x02, x08);
    r3  = L_add(x04, x06);
    r2  = L_sub(x04, x06);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y00 = L_add(x00, r1);
    r1  = L_add(y00, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    s1  = L_add(x03, x09);
    s4  = L_sub(x03, x09);
    s3  = L_add(x05, x07);
    s2  = L_sub(x05, x07);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y01 = L_add(x01, s1);
    s1  = L_add(y01, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y02 = L_add(r1, s2);
    y08 = L_sub(r1, s2);
    y04 = L_sub(r3, s4);
    y06 = L_add(r3, s4);

    y03 = L_sub(s1, r2);
    y09 = L_add(s1, r2);
    y05 = L_add(s3, r4);
    y07 = L_sub(s3, r4);

    /* 2. FFT5 stage */

    /* real part */
    r1  = L_add(x12, x18);
    r4  = L_sub(x12, x18);
    r3  = L_add(x14, x16);
    r2  = L_sub(x14, x16);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y10 = L_add(x10, r1);
    r1  = L_add(y10, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    s1  = L_add(x13, x19);
    s4  = L_sub(x13, x19);
    s3  = L_add(x15, x17);
    s2  = L_sub(x15, x17);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y11 = L_add(x11, s1);
    s1  = L_add(y11, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y12 = L_add(r1, s2);
    y18 = L_sub(r1, s2);
    y14 = L_sub(r3, s4);
    y16 = L_add(r3, s4);

    y13 = L_sub(s1, r2);
    y19 = L_add(s1, r2);
    y15 = L_add(s3, r4);
    y17 = L_sub(s3, r4);

    /* 3. FFT5 stage */

    /* real part */
    r1  = L_add(x22, x28);
    r4  = L_sub(x22, x28);
    r3  = L_add(x24, x26);
    r2  = L_sub(x24, x26);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y20 = L_add(x20, r1);
    r1  = L_add(y20, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    s1  = L_add(x23, x29);
    s4  = L_sub(x23, x29);
    s3  = L_add(x25, x27);
    s2  = L_sub(x25, x27);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y21 = L_add(x21, s1);
    s1  = L_add(y21, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y22 = L_add(r1, s2);
    y28 = L_sub(r1, s2);
    y24 = L_sub(r3, s4);
    y26 = L_add(r3, s4);

    y23 = L_sub(s1, r2);
    y29 = L_add(s1, r2);
    y25 = L_add(s3, r4);
    y27 = L_sub(s3, r4);

    /* 1. FFT3 stage */

    /* real part */
    r1  = L_add(y10, y20);
    r2  = Mpy_32_xx(L_sub(y10, y20), C31);
    z00 = L_add(y00, r1);
    r1  = L_sub(y00, L_shr_pos(r1, 1));

    /* imaginary part */
    s1  = L_add(y11, y21);
    s2  = Mpy_32_xx(L_sub(y11, y21), C31);
    z01 = L_add(y01, s1);
    s1  = L_sub(y01, L_shr_pos(s1, 1));

    /* combination */
    z20 = L_sub(r1, s2);
    z10 = L_add(r1, s2);
    z21 = L_add(s1, r2);
    z11 = L_sub(s1, r2);

    /* 2. FFT3 stage */

    /* real part */
    r1  = L_add(y12, y22);
    r2  = Mpy_32_xx(L_sub(y12, y22), C31);
    z12 = L_add(y02, r1);
    r1  = L_sub(y02, L_shr_pos(r1, 1));

    /* imaginary part */
    s1  = L_add(y13, y23);
    s2  = Mpy_32_xx(L_sub(y13, y23), C31);
    z13 = L_add(y03, s1);
    s1  = L_sub(y03, L_shr_pos(s1, 1));

    /* combination */
    z02 = L_sub(r1, s2);
    z22 = L_add(r1, s2);
    z03 = L_add(s1, r2);
    z23 = L_sub(s1, r2);

    /* 3. FFT3 stage */

    /* real part */
    r1  = L_add(y14, y24);
    r2  = Mpy_32_xx(L_sub(y14, y24), C31);
    z24 = L_add(y04, r1);
    r1  = L_sub(y04, L_shr_pos(r1, 1));

    /* imaginary part */
    s1  = L_add(y15, y25);
    s2  = Mpy_32_xx(L_sub(y15, y25), C31);
    z25 = L_add(y05, s1);
    s1  = L_sub(y05, L_shr_pos(s1, 1));

    /* combination */
    z14 = L_sub(r1, s2);
    z04 = L_add(r1, s2);
    z15 = L_add(s1, r2);
    z05 = L_sub(s1, r2);

    /* 4. FFT3 stage */

    /* real part */
    r1  = L_add(y16, y26);
    r2  = Mpy_32_xx(L_sub(y16, y26), C31);
    z06 = L_add(y06, r1);
    r1  = L_sub(y06, L_shr_pos(r1, 1));

    /* imaginary part */
    s1  = L_add(y17, y27);
    s2  = Mpy_32_xx(L_sub(y17, y27), C31);
    z07 = L_add(y07, s1);
    s1  = L_sub(y07, L_shr_pos(s1, 1));

    /* combination */
    z26 = L_sub(r1, s2);
    z16 = L_add(r1, s2);
    z27 = L_add(s1, r2);
    z17 = L_sub(s1, r2);

    /* 5. FFT3 stage */

    /* real part */
    r1  = L_add(y18, y28);
    r2  = Mpy_32_xx(L_sub(y18, y28), C31);
    z18 = L_add(y08, r1);
    r1  = L_sub(y08, L_shr_pos(r1, 1));

    /* imaginary part */
    s1  = L_add(y19, y29);
    s2  = Mpy_32_xx(L_sub(y19, y29), C31);
    z19 = L_add(y09, s1);
    s1  = L_sub(y09, L_shr_pos(s1, 1));

    /* combination */
    z08 = L_sub(r1, s2);
    z28 = L_add(r1, s2);
    z09 = L_add(s1, r2);
    z29 = L_sub(s1, r2);

    /* 2. FFT15 stage */
    x00 = L_shr_pos(re[s * 15], SCALEFACTOR30_1);
    x01 = L_shr_pos(im[s * 15], SCALEFACTOR30_1);
    x02 = L_shr_pos(re[s * 3], SCALEFACTOR30_1);
    x03 = L_shr_pos(im[s * 3], SCALEFACTOR30_1);
    x04 = L_shr_pos(re[s * 21], SCALEFACTOR30_1);
    x05 = L_shr_pos(im[s * 21], SCALEFACTOR30_1);
    x06 = L_shr_pos(re[s * 9], SCALEFACTOR30_1);
    x07 = L_shr_pos(im[s * 9], SCALEFACTOR30_1);
    x08 = L_shr_pos(re[s * 27], SCALEFACTOR30_1);
    x09 = L_shr_pos(im[s * 27], SCALEFACTOR30_1);

    x10 = L_shr_pos(re[s * 5], SCALEFACTOR30_1);
    x11 = L_shr_pos(im[s * 5], SCALEFACTOR30_1);
    x12 = L_shr_pos(re[s * 23], SCALEFACTOR30_1);
    x13 = L_shr_pos(im[s * 23], SCALEFACTOR30_1);
    x14 = L_shr_pos(re[s * 11], SCALEFACTOR30_1);
    x15 = L_shr_pos(im[s * 11], SCALEFACTOR30_1);
    x16 = L_shr_pos(re[s * 29], SCALEFACTOR30_1);
    x17 = L_shr_pos(im[s * 29], SCALEFACTOR30_1);
    x18 = L_shr_pos(re[s * 17], SCALEFACTOR30_1);
    x19 = L_shr_pos(im[s * 17], SCALEFACTOR30_1);

    x20 = L_shr_pos(re[s * 25], SCALEFACTOR30_1);
    x21 = L_shr_pos(im[s * 25], SCALEFACTOR30_1);
    x22 = L_shr_pos(re[s * 13], SCALEFACTOR30_1);
    x23 = L_shr_pos(im[s * 13], SCALEFACTOR30_1);
    x24 = L_shr_pos(re[s * 1], SCALEFACTOR30_1);
    x25 = L_shr_pos(im[s * 1], SCALEFACTOR30_1);
    x26 = L_shr_pos(re[s * 19], SCALEFACTOR30_1);
    x27 = L_shr_pos(im[s * 19], SCALEFACTOR30_1);
    x28 = L_shr_pos(re[s * 7], SCALEFACTOR30_1);
    x29 = L_shr_pos(im[s * 7], SCALEFACTOR30_1);

    /* 1. FFT5 stage */

    /* real part */
    r1  = L_add(x02, x08);
    r4  = L_sub(x02, x08);
    r3  = L_add(x04, x06);
    r2  = L_sub(x04, x06);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y00 = L_add(x00, r1);
    r1  = L_add(y00, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    s1  = L_add(x03, x09);
    s4  = L_sub(x03, x09);
    s3  = L_add(x05, x07);
    s2  = L_sub(x05, x07);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y01 = L_add(x01, s1);
    s1  = L_add(y01, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y02 = L_add(r1, s2);
    y08 = L_sub(r1, s2);
    y04 = L_sub(r3, s4);
    y06 = L_add(r3, s4);

    y03 = L_sub(s1, r2);
    y09 = L_add(s1, r2);
    y05 = L_add(s3, r4);
    y07 = L_sub(s3, r4);

    /* 2. FFT5 stage */

    /* real part */
    r1  = L_add(x12, x18);
    r4  = L_sub(x12, x18);
    r3  = L_add(x14, x16);
    r2  = L_sub(x14, x16);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y10 = L_add(x10, r1);
    r1  = L_add(y10, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    s1  = L_add(x13, x19);
    s4  = L_sub(x13, x19);
    s3  = L_add(x15, x17);
    s2  = L_sub(x15, x17);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y11 = L_add(x11, s1);
    s1  = L_add(y11, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y12 = L_add(r1, s2);
    y18 = L_sub(r1, s2);
    y14 = L_sub(r3, s4);
    y16 = L_add(r3, s4);

    y13 = L_sub(s1, r2);
    y19 = L_add(s1, r2);
    y15 = L_add(s3, r4);
    y17 = L_sub(s3, r4);

    /* 3. FFT5 stage */

    /* real part */
    r1  = L_add(x22, x28);
    r4  = L_sub(x22, x28);
    r3  = L_add(x24, x26);
    r2  = L_sub(x24, x26);
    t   = Mpy_32_xx(L_sub(r1, r3), C54);
    r1  = L_add(r1, r3);
    y20 = L_add(x20, r1);
    r1  = L_add(y20, (L_shl_pos(Mpy_32_xx(r1, C55), 1)));
    r3  = L_sub(r1, t);
    r1  = L_add(r1, t);
    t   = Mpy_32_xx((L_add(r4, r2)), C51);
    r4  = L_add(t, L_shl_pos(Mpy_32_xx(r4, C52), 1));
    r2  = L_add(t, Mpy_32_xx(r2, C53));

    /* imaginary part */
    s1  = L_add(x23, x29);
    s4  = L_sub(x23, x29);
    s3  = L_add(x25, x27);
    s2  = L_sub(x25, x27);
    t   = Mpy_32_xx(L_sub(s1, s3), C54);
    s1  = L_add(s1, s3);
    y21 = L_add(x21, s1);
    s1  = L_add(y21, L_shl_pos(Mpy_32_xx(s1, C55), 1));
    s3  = L_sub(s1, t);
    s1  = L_add(s1, t);
    t   = Mpy_32_xx(L_add(s4, s2), C51);
    s4  = L_add(t, L_shl_pos(Mpy_32_xx(s4, C52), 1));
    s2  = L_add(t, Mpy_32_xx(s2, C53));

    /* combination */
    y22 = L_add(r1, s2);
    y28 = L_sub(r1, s2);
    y24 = L_sub(r3, s4);
    y26 = L_add(r3, s4);

    y23 = L_sub(s1, r2);
    y29 = L_add(s1, r2);
    y25 = L_add(s3, r4);
    y27 = L_sub(s3, r4);

    /* 1. FFT3 stage */

    /* real part */
    r1  = L_add(y10, y20);
    r2  = Mpy_32_xx(L_sub(y10, y20), C31);
    z30 = L_add(y00, r1);
    r1  = L_sub(y00, L_shr_pos(r1, 1));

    /* imaginary part */
    s1  = L_add(y11, y21);
    s2  = Mpy_32_xx(L_sub(y11, y21), C31);
    z31 = L_add(y01, s1);
    s1  = L_sub(y01, L_shr_pos(s1, 1));

    /* combination */
    z50 = L_sub(r1, s2);
    z40 = L_add(r1, s2);
    z51 = L_add(s1, r2);
    z41 = L_sub(s1, r2);

    /* 2. FFT3 stage */

    /* real part */
    r1  = L_add(y12, y22);
    r2  = Mpy_32_xx(L_sub(y12, y22), C31);
    z42 = L_add(y02, r1);
    r1  = L_sub(y02, L_shr_pos(r1, 1));

    /* imaginary part */
    s1  = L_add(y13, y23);
    s2  = Mpy_32_xx(L_sub(y13, y23), C31);
    z43 = L_add(y03, s1);
    s1  = L_sub(y03, L_shr_pos(s1, 1));

    /* combination */
    z32 = L_sub(r1, s2);
    z52 = L_add(r1, s2);
    z33 = L_add(s1, r2);
    z53 = L_sub(s1, r2);

    /* 3. FFT3 stage */

    /* real part */
    r1  = L_add(y14, y24);
    r2  = Mpy_32_xx(L_sub(y14, y24), C31);
    z54 = L_add(y04, r1);
    r1  = L_sub(y04, L_shr_pos(r1, 1));

    /* imaginary part */
    s1  = L_add(y15, y25);
    s2  = Mpy_32_xx(L_sub(y15, y25), C31);
    z55 = L_add(y05, s1);
    s1  = L_sub(y05, L_shr_pos(s1, 1));

    /* combination */
    z44 = L_sub(r1, s2);
    z34 = L_add(r1, s2);
    z45 = L_add(s1, r2);
    z35 = L_sub(s1, r2);

    /* 4. FFT3 stage */

    /* real part */
    r1  = L_add(y16, y26);
    r2  = Mpy_32_xx(L_sub(y16, y26), C31);
    z36 = L_add(y06, r1);
    r1  = L_sub(y06, L_shr_pos(r1, 1));

    /* imaginary part */
    s1  = L_add(y17, y27);
    s2  = Mpy_32_xx(L_sub(y17, y27), C31);
    z37 = L_add(y07, s1);
    s1  = L_sub(y07, L_shr_pos(s1, 1));

    /* combination */
    z56 = L_sub(r1, s2);
    z46 = L_add(r1, s2);
    z57 = L_add(s1, r2);
    z47 = L_sub(s1, r2);

    /* 5. FFT3 stage */

    /* real part */
    r1  = L_add(y18, y28);
    r2  = Mpy_32_xx(L_sub(y18, y28), C31);
    z48 = L_add(y08, r1);
    r1  = L_sub(y08, L_shr_pos(r1, 1));

    /* imaginary part */
    s1  = L_add(y19, y29);
    s2  = Mpy_32_xx(L_sub(y19, y29), C31);
    z49 = L_add(y09, s1);
    s1  = L_sub(y09, L_shr_pos(s1, 1));

    /* combination */
    z38 = L_sub(r1, s2);
    z58 = L_add(r1, s2);
    z39 = L_add(s1, r2);
    z59 = L_sub(s1, r2);

    /* 1. FFT2 stage */
    r1   = L_shr_pos(z00, SCALEFACTOR30_2);
    r2   = L_shr_pos(z30, SCALEFACTOR30_2);
    r3   = L_shr_pos(z01, SCALEFACTOR30_2);
    r4   = L_shr_pos(z31, SCALEFACTOR30_2);
    *rel = L_add(r1, r2); move32();
    *reh = L_sub(r1, r2); move32();
    *iml = L_add(r3, r4); move32();
    *imh = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 2. FFT2 stage */
    r1   = L_shr_pos(z16, SCALEFACTOR30_2);
    r2   = L_shr_pos(z46, SCALEFACTOR30_2);
    r3   = L_shr_pos(z17, SCALEFACTOR30_2);
    r4   = L_shr_pos(z47, SCALEFACTOR30_2);
    *reh = L_add(r1, r2); move32();
    *rel = L_sub(r1, r2); move32();
    *imh = L_add(r3, r4); move32();
    *iml = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 3. FFT2 stage */
    r1   = L_shr_pos(z02, SCALEFACTOR30_2);
    r2   = L_shr_pos(z32, SCALEFACTOR30_2);
    r3   = L_shr_pos(z03, SCALEFACTOR30_2);
    r4   = L_shr_pos(z33, SCALEFACTOR30_2);
    *rel = L_add(r1, r2); move32();
    *reh = L_sub(r1, r2); move32();
    *iml = L_add(r3, r4); move32();
    *imh = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 4. FFT2 stage */
    r1   = L_shr_pos(z18, SCALEFACTOR30_2);
    r2   = L_shr_pos(z48, SCALEFACTOR30_2);
    r3   = L_shr_pos(z19, SCALEFACTOR30_2);
    r4   = L_shr_pos(z49, SCALEFACTOR30_2);
    *reh = L_add(r1, r2); move32();
    *rel = L_sub(r1, r2); move32();
    *imh = L_add(r3, r4); move32();
    *iml = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 5. FFT2 stage */
    r1   = L_shr_pos(z04, SCALEFACTOR30_2);
    r2   = L_shr_pos(z34, SCALEFACTOR30_2);
    r3   = L_shr_pos(z05, SCALEFACTOR30_2);
    r4   = L_shr_pos(z35, SCALEFACTOR30_2);
    *rel = L_add(r1, r2); move32();
    *reh = L_sub(r1, r2); move32();
    *iml = L_add(r3, r4); move32();
    *imh = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 6. FFT2 stage */
    r1   = L_shr_pos(z20, SCALEFACTOR30_2);
    r2   = L_shr_pos(z50, SCALEFACTOR30_2);
    r3   = L_shr_pos(z21, SCALEFACTOR30_2);
    r4   = L_shr_pos(z51, SCALEFACTOR30_2);
    *reh = L_add(r1, r2); move32();
    *rel = L_sub(r1, r2); move32();
    *imh = L_add(r3, r4); move32();
    *iml = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 7. FFT2 stage */
    r1   = L_shr_pos(z06, SCALEFACTOR30_2);
    r2   = L_shr_pos(z36, SCALEFACTOR30_2);
    r3   = L_shr_pos(z07, SCALEFACTOR30_2);
    r4   = L_shr_pos(z37, SCALEFACTOR30_2);
    *rel = L_add(r1, r2); move32();
    *reh = L_sub(r1, r2); move32();
    *iml = L_add(r3, r4); move32();
    *imh = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 8. FFT2 stage */
    r1   = L_shr_pos(z22, SCALEFACTOR30_2);
    r2   = L_shr_pos(z52, SCALEFACTOR30_2);
    r3   = L_shr_pos(z23, SCALEFACTOR30_2);
    r4   = L_shr_pos(z53, SCALEFACTOR30_2);
    *reh = L_add(r1, r2); move32();
    *rel = L_sub(r1, r2); move32();
    *imh = L_add(r3, r4); move32();
    *iml = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 9. FFT2 stage */
    r1   = L_shr_pos(z08, SCALEFACTOR30_2);
    r2   = L_shr_pos(z38, SCALEFACTOR30_2);
    r3   = L_shr_pos(z09, SCALEFACTOR30_2);
    r4   = L_shr_pos(z39, SCALEFACTOR30_2);
    *rel = L_add(r1, r2); move32();
    *reh = L_sub(r1, r2); move32();
    *iml = L_add(r3, r4); move32();
    *imh = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 10. FFT2 stage */
    r1   = L_shr_pos(z24, SCALEFACTOR30_2);
    r2   = L_shr_pos(z54, SCALEFACTOR30_2);
    r3   = L_shr_pos(z25, SCALEFACTOR30_2);
    r4   = L_shr_pos(z55, SCALEFACTOR30_2);
    *reh = L_add(r1, r2); move32();
    *rel = L_sub(r1, r2); move32();
    *imh = L_add(r3, r4); move32();
    *iml = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 11. FFT2 stage */
    r1   = L_shr_pos(z10, SCALEFACTOR30_2);
    r2   = L_shr_pos(z40, SCALEFACTOR30_2);
    r3   = L_shr_pos(z11, SCALEFACTOR30_2);
    r4   = L_shr_pos(z41, SCALEFACTOR30_2);
    *rel = L_add(r1, r2); move32();
    *reh = L_sub(r1, r2); move32();
    *iml = L_add(r3, r4); move32();
    *imh = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 12. FFT2 stage */
    r1   = L_shr_pos(z26, SCALEFACTOR30_2);
    r2   = L_shr_pos(z56, SCALEFACTOR30_2);
    r3   = L_shr_pos(z27, SCALEFACTOR30_2);
    r4   = L_shr_pos(z57, SCALEFACTOR30_2);
    *reh = L_add(r1, r2); move32();
    *rel = L_sub(r1, r2); move32();
    *imh = L_add(r3, r4); move32();
    *iml = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 13. FFT2 stage */
    r1   = L_shr_pos(z12, SCALEFACTOR30_2);
    r2   = L_shr_pos(z42, SCALEFACTOR30_2);
    r3   = L_shr_pos(z13, SCALEFACTOR30_2);
    r4   = L_shr_pos(z43, SCALEFACTOR30_2);
    *rel = L_add(r1, r2); move32();
    *reh = L_sub(r1, r2); move32();
    *iml = L_add(r3, r4); move32();
    *imh = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 14. FFT2 stage */
    r1   = L_shr_pos(z28, SCALEFACTOR30_2);
    r2   = L_shr_pos(z58, SCALEFACTOR30_2);
    r3   = L_shr_pos(z29, SCALEFACTOR30_2);
    r4   = L_shr_pos(z59, SCALEFACTOR30_2);
    *reh = L_add(r1, r2); move32();
    *rel = L_sub(r1, r2); move32();
    *imh = L_add(r3, r4); move32();
    *iml = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    /* 15. FFT2 stage */
    r1   = L_shr_pos(z14, SCALEFACTOR30_2);
    r2   = L_shr_pos(z44, SCALEFACTOR30_2);
    r3   = L_shr_pos(z15, SCALEFACTOR30_2);
    r4   = L_shr_pos(z45, SCALEFACTOR30_2);
    *rel = L_add(r1, r2); move32();
    *reh = L_sub(r1, r2); move32();
    *iml = L_add(r3, r4); move32();
    *imh = L_sub(r3, r4); move32();
    rel += s, reh += s, iml += s;
    imh += s;

    Dyn_Mem_Deluxe_Out();
}

/**
 * \brief    Function performs a complex 32-point FFT
 *           The FFT is performed inplace. The result of the FFT
 *           is scaled by SCALEFACTOR32 bits.
 *
 *           WOPS with 32x16 bit multiplications:  752 cycles
 *
 * \param    [i/o] re    real input / output
 * \param    [i/o] im    imag input / output
 * \param    [i  ] s     stride real and imag input / output
 *
 * \return   void
 */


static void fft32(Word32 *re, Word32 *im, Word16 s)
{
    Dyn_Mem_Deluxe_In(
        Word32 as, bs;
        Word32 x00, x01, x02, x03, x04, x05, x06, x07;
        Word32 x08, x09, x10, x11, x12, x13, x14, x15;
        Word32 t00, t01, t02, t03, t04, t05, t06, t07;
        Word32 t08, t09, t10, t11, t12, t13, t14, t15;
        Word32 s00, s01, s02, s03, s04, s05, s06, s07;
        Word32 s08, s09, s10, s11, s12, s13, s14, s15;

        Word32 y00, y01, y02, y03, y04, y05, y06, y07;
        Word32 y08, y09, y10, y11, y12, y13, y14, y15;
        Word32 y16, y17, y18, y19, y20, y21, y22, y23;
        Word32 y24, y25, y26, y27, y28, y29, y30, y31;
        Word32 y32, y33, y34, y35, y36, y37, y38, y39;
        Word32 y40, y41, y42, y43, y44, y45, y46, y47;
        Word32 y48, y49, y50, y51, y52, y53, y54, y55;
        Word32 y56, y57, y58, y59, y60, y61, y62, y63;
    );

    /* 1. FFT8 stage */
    x00 = L_shr_pos(re[s * 0], SCALEFACTOR32_1);
    x01 = L_shr_pos(im[s * 0], SCALEFACTOR32_1);
    x02 = L_shr_pos(re[s * 4], SCALEFACTOR32_1);
    x03 = L_shr_pos(im[s * 4], SCALEFACTOR32_1);
    x04 = L_shr_pos(re[s * 8], SCALEFACTOR32_1);
    x05 = L_shr_pos(im[s * 8], SCALEFACTOR32_1);
    x06 = L_shr_pos(re[s * 12], SCALEFACTOR32_1);
    x07 = L_shr_pos(im[s * 12], SCALEFACTOR32_1);
    x08 = L_shr_pos(re[s * 16], SCALEFACTOR32_1);
    x09 = L_shr_pos(im[s * 16], SCALEFACTOR32_1);
    x10 = L_shr_pos(re[s * 20], SCALEFACTOR32_1);
    x11 = L_shr_pos(im[s * 20], SCALEFACTOR32_1);
    x12 = L_shr_pos(re[s * 24], SCALEFACTOR32_1);
    x13 = L_shr_pos(im[s * 24], SCALEFACTOR32_1);
    x14 = L_shr_pos(re[s * 28], SCALEFACTOR32_1);
    x15 = L_shr_pos(im[s * 28], SCALEFACTOR32_1);

    t00 = L_add(x00, x08);
    t02 = L_sub(x00, x08);
    t01 = L_add(x01, x09);
    t03 = L_sub(x01, x09);
    t04 = L_add(x02, x10);
    t06 = L_sub(x02, x10);
    t05 = L_add(x03, x11);
    t07 = L_sub(x03, x11);
    t08 = L_add(x04, x12);
    t10 = L_sub(x04, x12);
    t09 = L_add(x05, x13);
    t11 = L_sub(x05, x13);
    t12 = L_add(x06, x14);
    t14 = L_sub(x06, x14);
    t13 = L_add(x07, x15);
    t15 = L_sub(x07, x15);

    /* Pre-additions and core multiplications */
    s00 = L_add(t00, t08);
    s04 = L_sub(t00, t08);
    s01 = L_add(t01, t09);
    s05 = L_sub(t01, t09);
    s08 = L_sub(t02, t11);
    s10 = L_add(t02, t11);
    s09 = L_add(t03, t10);
    s11 = L_sub(t03, t10);
    s02 = L_add(t04, t12);
    s07 = L_sub(t04, t12);
    s03 = L_add(t05, t13);
    s06 = L_sub(t13, t05);
    t01 = L_add(t06, t14);
    t02 = L_sub(t06, t14);
    t00 = L_add(t07, t15);
    t03 = L_sub(t07, t15);

    Mpy3_0(s12, s13, s14, s15, t00, t01, t02, t03);

    /* Post-additions */
    y00 = L_add(s00, s02);
    y08 = L_sub(s00, s02);
    y01 = L_add(s01, s03);
    y09 = L_sub(s01, s03);
    y04 = L_sub(s04, s06);
    y12 = L_add(s04, s06);
    y05 = L_sub(s05, s07);
    y13 = L_add(s05, s07);
    y06 = L_add(s08, s14);
    y14 = L_sub(s08, s14);
    y07 = L_add(s09, s15);
    y15 = L_sub(s09, s15);
    y02 = L_add(s10, s12);
    y10 = L_sub(s10, s12);
    y03 = L_add(s11, s13);
    y11 = L_sub(s11, s13);

    /* 2. FFT8 stage */
    x00 = L_shr_pos(re[s * 1], SCALEFACTOR32_1);
    x01 = L_shr_pos(im[s * 1], SCALEFACTOR32_1);
    x02 = L_shr_pos(re[s * 5], SCALEFACTOR32_1);
    x03 = L_shr_pos(im[s * 5], SCALEFACTOR32_1);
    x04 = L_shr_pos(re[s * 9], SCALEFACTOR32_1);
    x05 = L_shr_pos(im[s * 9], SCALEFACTOR32_1);
    x06 = L_shr_pos(re[s * 13], SCALEFACTOR32_1);
    x07 = L_shr_pos(im[s * 13], SCALEFACTOR32_1);
    x08 = L_shr_pos(re[s * 17], SCALEFACTOR32_1);
    x09 = L_shr_pos(im[s * 17], SCALEFACTOR32_1);
    x10 = L_shr_pos(re[s * 21], SCALEFACTOR32_1);
    x11 = L_shr_pos(im[s * 21], SCALEFACTOR32_1);
    x12 = L_shr_pos(re[s * 25], SCALEFACTOR32_1);
    x13 = L_shr_pos(im[s * 25], SCALEFACTOR32_1);
    x14 = L_shr_pos(re[s * 29], SCALEFACTOR32_1);
    x15 = L_shr_pos(im[s * 29], SCALEFACTOR32_1);

    t00 = L_add(x00, x08);
    t02 = L_sub(x00, x08);
    t01 = L_add(x01, x09);
    t03 = L_sub(x01, x09);
    t04 = L_add(x02, x10);
    t06 = L_sub(x02, x10);
    t05 = L_add(x03, x11);
    t07 = L_sub(x03, x11);
    t08 = L_add(x04, x12);
    t10 = L_sub(x04, x12);
    t09 = L_add(x05, x13);
    t11 = L_sub(x05, x13);
    t12 = L_add(x06, x14);
    t14 = L_sub(x06, x14);
    t13 = L_add(x07, x15);
    t15 = L_sub(x07, x15);

    /* Pre-additions and core multiplications */
    s00 = L_add(t00, t08);
    s04 = L_sub(t00, t08);
    s01 = L_add(t01, t09);
    s05 = L_sub(t01, t09);
    s08 = L_sub(t02, t11);
    s10 = L_add(t02, t11);
    s09 = L_add(t03, t10);
    s11 = L_sub(t03, t10);
    s02 = L_add(t04, t12);
    s07 = L_sub(t04, t12);
    s03 = L_add(t05, t13);
    s06 = L_sub(t13, t05);
    t01 = L_add(t06, t14);
    t02 = L_sub(t06, t14);
    t00 = L_add(t07, t15);
    t03 = L_sub(t07, t15);

    Mpy3_0(s12, s13, s14, s15, t00, t01, t02, t03);

    /* Post-additions */
    y16 = L_add(s00, s02);
    y24 = L_sub(s00, s02);
    y17 = L_add(s01, s03);
    y25 = L_sub(s01, s03);
    y20 = L_sub(s04, s06);
    y28 = L_add(s04, s06);
    y21 = L_sub(s05, s07);
    y29 = L_add(s05, s07);
    y22 = L_add(s08, s14);
    y30 = L_sub(s08, s14);
    y23 = L_add(s09, s15);
    y31 = L_sub(s09, s15);
    y18 = L_add(s10, s12);
    y26 = L_sub(s10, s12);
    y19 = L_add(s11, s13);
    y27 = L_sub(s11, s13);

    /* 3. FFT8 stage */
    x00 = L_shr_pos(re[s * 2], SCALEFACTOR32_1);
    x01 = L_shr_pos(im[s * 2], SCALEFACTOR32_1);
    x02 = L_shr_pos(re[s * 6], SCALEFACTOR32_1);
    x03 = L_shr_pos(im[s * 6], SCALEFACTOR32_1);
    x04 = L_shr_pos(re[s * 10], SCALEFACTOR32_1);
    x05 = L_shr_pos(im[s * 10], SCALEFACTOR32_1);
    x06 = L_shr_pos(re[s * 14], SCALEFACTOR32_1);
    x07 = L_shr_pos(im[s * 14], SCALEFACTOR32_1);
    x08 = L_shr_pos(re[s * 18], SCALEFACTOR32_1);
    x09 = L_shr_pos(im[s * 18], SCALEFACTOR32_1);
    x10 = L_shr_pos(re[s * 22], SCALEFACTOR32_1);
    x11 = L_shr_pos(im[s * 22], SCALEFACTOR32_1);
    x12 = L_shr_pos(re[s * 26], SCALEFACTOR32_1);
    x13 = L_shr_pos(im[s * 26], SCALEFACTOR32_1);
    x14 = L_shr_pos(re[s * 30], SCALEFACTOR32_1);
    x15 = L_shr_pos(im[s * 30], SCALEFACTOR32_1);

    t00 = L_add(x00, x08);
    t02 = L_sub(x00, x08);
    t01 = L_add(x01, x09);
    t03 = L_sub(x01, x09);
    t04 = L_add(x02, x10);
    t06 = L_sub(x02, x10);
    t05 = L_add(x03, x11);
    t07 = L_sub(x03, x11);
    t08 = L_add(x04, x12);
    t10 = L_sub(x04, x12);
    t09 = L_add(x05, x13);
    t11 = L_sub(x05, x13);
    t12 = L_add(x06, x14);
    t14 = L_sub(x06, x14);
    t13 = L_add(x07, x15);
    t15 = L_sub(x07, x15);

    /* Pre-additions and core multiplications */
    s00 = L_add(t00, t08);
    s04 = L_sub(t00, t08);
    s01 = L_add(t01, t09);
    s05 = L_sub(t01, t09);
    s08 = L_sub(t02, t11);
    s10 = L_add(t02, t11);
    s09 = L_add(t03, t10);
    s11 = L_sub(t03, t10);
    s02 = L_add(t04, t12);
    s07 = L_sub(t04, t12);
    s03 = L_add(t05, t13);
    s06 = L_sub(t13, t05);
    t01 = L_add(t06, t14);
    t02 = L_sub(t06, t14);
    t00 = L_add(t07, t15);
    t03 = L_sub(t07, t15);

    Mpy3_0(s12, s13, s14, s15, t00, t01, t02, t03);

    /* Post-additions */
    y32 = L_add(s00, s02);
    y40 = L_sub(s00, s02);
    y33 = L_add(s01, s03);
    y41 = L_sub(s01, s03);
    y36 = L_sub(s04, s06);
    y44 = L_add(s04, s06);
    y37 = L_sub(s05, s07);
    y45 = L_add(s05, s07);
    y38 = L_add(s08, s14);
    y46 = L_sub(s08, s14);
    y39 = L_add(s09, s15);
    y47 = L_sub(s09, s15);
    y34 = L_add(s10, s12);
    y42 = L_sub(s10, s12);
    y35 = L_add(s11, s13);
    y43 = L_sub(s11, s13);

    /* 4. FFT8 stage */
    x00 = L_shr_pos(re[s * 3], SCALEFACTOR32_1);
    x01 = L_shr_pos(im[s * 3], SCALEFACTOR32_1);
    x02 = L_shr_pos(re[s * 7], SCALEFACTOR32_1);
    x03 = L_shr_pos(im[s * 7], SCALEFACTOR32_1);
    x04 = L_shr_pos(re[s * 11], SCALEFACTOR32_1);
    x05 = L_shr_pos(im[s * 11], SCALEFACTOR32_1);
    x06 = L_shr_pos(re[s * 15], SCALEFACTOR32_1);
    x07 = L_shr_pos(im[s * 15], SCALEFACTOR32_1);
    x08 = L_shr_pos(re[s * 19], SCALEFACTOR32_1);
    x09 = L_shr_pos(im[s * 19], SCALEFACTOR32_1);
    x10 = L_shr_pos(re[s * 23], SCALEFACTOR32_1);
    x11 = L_shr_pos(im[s * 23], SCALEFACTOR32_1);
    x12 = L_shr_pos(re[s * 27], SCALEFACTOR32_1);
    x13 = L_shr_pos(im[s * 27], SCALEFACTOR32_1);
    x14 = L_shr_pos(re[s * 31], SCALEFACTOR32_1);
    x15 = L_shr_pos(im[s * 31], SCALEFACTOR32_1);

    t00 = L_add(x00, x08);
    t02 = L_sub(x00, x08);
    t01 = L_add(x01, x09);
    t03 = L_sub(x01, x09);
    t04 = L_add(x02, x10);
    t06 = L_sub(x02, x10);
    t05 = L_add(x03, x11);
    t07 = L_sub(x03, x11);
    t08 = L_add(x04, x12);
    t10 = L_sub(x04, x12);
    t09 = L_add(x05, x13);
    t11 = L_sub(x05, x13);
    t12 = L_add(x06, x14);
    t14 = L_sub(x06, x14);
    t13 = L_add(x07, x15);
    t15 = L_sub(x07, x15);

    /* Pre-additions and core multiplications */
    s00 = L_add(t00, t08);
    s04 = L_sub(t00, t08);
    s01 = L_add(t01, t09);
    s05 = L_sub(t01, t09);
    s08 = L_sub(t02, t11);
    s10 = L_add(t02, t11);
    s09 = L_add(t03, t10);
    s11 = L_sub(t03, t10);
    s02 = L_add(t04, t12);
    s07 = L_sub(t04, t12);
    s03 = L_add(t05, t13);
    s06 = L_sub(t13, t05);
    t01 = L_add(t06, t14);
    t02 = L_sub(t06, t14);
    t00 = L_add(t07, t15);
    t03 = L_sub(t07, t15);

    Mpy3_0(s12, s13, s14, s15, t00, t01, t02, t03);

    /* Post-additions */
    y48 = L_add(s00, s02);
    y56 = L_sub(s00, s02);
    y49 = L_add(s01, s03);
    y57 = L_sub(s01, s03);
    y52 = L_sub(s04, s06);
    y60 = L_add(s04, s06);
    y53 = L_sub(s05, s07);
    y61 = L_add(s05, s07);
    y54 = L_add(s08, s14);
    y62 = L_sub(s08, s14);
    y55 = L_add(s09, s15);
    y63 = L_sub(s09, s15);
    y50 = L_add(s10, s12);
    y58 = L_sub(s10, s12);
    y51 = L_add(s11, s13);
    y59 = L_sub(s11, s13);

    /* apply twiddle factors */
    y00 = L_shr_pos(y00, SCALEFACTOR32_2);
    y01 = L_shr_pos(y01, SCALEFACTOR32_2);
    y02 = L_shr_pos(y02, SCALEFACTOR32_2);
    y03 = L_shr_pos(y03, SCALEFACTOR32_2);
    y04 = L_shr_pos(y04, SCALEFACTOR32_2);
    y05 = L_shr_pos(y05, SCALEFACTOR32_2);
    y06 = L_shr_pos(y06, SCALEFACTOR32_2);
    y07 = L_shr_pos(y07, SCALEFACTOR32_2);
    y08 = L_shr_pos(y08, SCALEFACTOR32_2);
    y09 = L_shr_pos(y09, SCALEFACTOR32_2);
    y10 = L_shr_pos(y10, SCALEFACTOR32_2);
    y11 = L_shr_pos(y11, SCALEFACTOR32_2);
    y12 = L_shr_pos(y12, SCALEFACTOR32_2);
    y13 = L_shr_pos(y13, SCALEFACTOR32_2);
    y14 = L_shr_pos(y14, SCALEFACTOR32_2);
    y15 = L_shr_pos(y15, SCALEFACTOR32_2);
    y16 = L_shr_pos(y16, SCALEFACTOR32_2);
    y17 = L_shr_pos(y17, SCALEFACTOR32_2);
    y32 = L_shr_pos(y32, SCALEFACTOR32_2);
    y33 = L_shr_pos(y33, SCALEFACTOR32_2);
    y48 = L_shr_pos(y48, SCALEFACTOR32_2);
    y49 = L_shr_pos(y49, SCALEFACTOR32_2);
    y40 = L_shr_pos(y40, SCALEFACTOR32_2);
    y41 = L_shr_pos(y41, SCALEFACTOR32_2);

    cplxMpy3_0(y18, y19, RotVector_32_32[2 * 0 + 0], RotVector_32_32[2 * 0 + 1]);
    cplxMpy3_0(y20, y21, RotVector_32_32[2 * 1 + 0], RotVector_32_32[2 * 1 + 1]);
    cplxMpy3_0(y22, y23, RotVector_32_32[2 * 2 + 0], RotVector_32_32[2 * 2 + 1]);
    cplxMpy3_0(y24, y25, RotVector_32_32[2 * 3 + 0], RotVector_32_32[2 * 3 + 1]);
    cplxMpy3_0(y26, y27, RotVector_32_32[2 * 4 + 0], RotVector_32_32[2 * 4 + 1]);
    cplxMpy3_0(y28, y29, RotVector_32_32[2 * 5 + 0], RotVector_32_32[2 * 5 + 1]);
    cplxMpy3_0(y30, y31, RotVector_32_32[2 * 6 + 0], RotVector_32_32[2 * 6 + 1]);
    cplxMpy3_0(y34, y35, RotVector_32_32[2 * 7 + 0], RotVector_32_32[2 * 7 + 1]);
    cplxMpy3_0(y36, y37, RotVector_32_32[2 * 8 + 0], RotVector_32_32[2 * 8 + 1]);
    cplxMpy3_0(y38, y39, RotVector_32_32[2 * 9 + 0], RotVector_32_32[2 * 9 + 1]);
    cplxMpy3_0(y42, y43, RotVector_32_32[2 * 10 + 0], RotVector_32_32[2 * 10 + 1]);
    cplxMpy3_0(y44, y45, RotVector_32_32[2 * 11 + 0], RotVector_32_32[2 * 11 + 1]);
    cplxMpy3_0(y46, y47, RotVector_32_32[2 * 12 + 0], RotVector_32_32[2 * 12 + 1]);
    cplxMpy3_0(y50, y51, RotVector_32_32[2 * 13 + 0], RotVector_32_32[2 * 13 + 1]);
    cplxMpy3_0(y52, y53, RotVector_32_32[2 * 14 + 0], RotVector_32_32[2 * 14 + 1]);
    cplxMpy3_0(y54, y55, RotVector_32_32[2 * 15 + 0], RotVector_32_32[2 * 15 + 1]);
    cplxMpy3_0(y56, y57, RotVector_32_32[2 * 16 + 0], RotVector_32_32[2 * 16 + 1]);
    cplxMpy3_0(y58, y59, RotVector_32_32[2 * 17 + 0], RotVector_32_32[2 * 17 + 1]);
    cplxMpy3_0(y60, y61, RotVector_32_32[2 * 18 + 0], RotVector_32_32[2 * 18 + 1]);
    cplxMpy3_0(y62, y63, RotVector_32_32[2 * 19 + 0], RotVector_32_32[2 * 19 + 1]);

    /* 1. FFT4 stage */

    /* Pre-additions */
    t00 = L_add(y00, y32);
    t02 = L_sub(y00, y32);
    t01 = L_add(y01, y33);
    t03 = L_sub(y01, y33);
    t04 = L_add(y16, y48);
    t07 = L_sub(y16, y48);
    t05 = L_add(y49, y17);
    t06 = L_sub(y49, y17);

    /* Post-additions */
    re[s * 0]  = L_add(t00, t04); move32();
    im[s * 0]  = L_add(t01, t05); move32();
    re[s * 8]  = L_sub(t02, t06); move32();
    im[s * 8]  = L_sub(t03, t07); move32();
    re[s * 16] = L_sub(t00, t04); move32();
    im[s * 16] = L_sub(t01, t05); move32();
    re[s * 24] = L_add(t02, t06); move32();
    im[s * 24] = L_add(t03, t07); move32();

    /* 2. FFT4 stage */

    /* Pre-additions */
    t00 = L_add(y02, y34);
    t02 = L_sub(y02, y34);
    t01 = L_add(y03, y35);
    t03 = L_sub(y03, y35);
    t04 = L_add(y18, y50);
    t07 = L_sub(y18, y50);
    t05 = L_add(y51, y19);
    t06 = L_sub(y51, y19);

    /* Post-additions */
    re[s * 1]  = L_add(t00, t04); move32();
    im[s * 1]  = L_add(t01, t05); move32();
    re[s * 9]  = L_sub(t02, t06); move32();
    im[s * 9]  = L_sub(t03, t07); move32();
    re[s * 17] = L_sub(t00, t04); move32();
    im[s * 17] = L_sub(t01, t05); move32();
    re[s * 25] = L_add(t02, t06); move32();
    im[s * 25] = L_add(t03, t07); move32();

    /* 3. FFT4 stage */

    /* Pre-additions */
    t00 = L_add(y04, y36);
    t02 = L_sub(y04, y36);
    t01 = L_add(y05, y37);
    t03 = L_sub(y05, y37);
    t04 = L_add(y20, y52);
    t07 = L_sub(y20, y52);
    t05 = L_add(y53, y21);
    t06 = L_sub(y53, y21);

    /* Post-additions */
    re[s * 2]  = L_add(t00, t04); move32();
    im[s * 2]  = L_add(t01, t05); move32();
    re[s * 10] = L_sub(t02, t06); move32();
    im[s * 10] = L_sub(t03, t07); move32();
    re[s * 18] = L_sub(t00, t04); move32();
    im[s * 18] = L_sub(t01, t05); move32();
    re[s * 26] = L_add(t02, t06); move32();
    im[s * 26] = L_add(t03, t07); move32();

    /* 4. FFT4 stage */

    /* Pre-additions */
    t00 = L_add(y06, y38);
    t02 = L_sub(y06, y38);
    t01 = L_add(y07, y39);
    t03 = L_sub(y07, y39);
    t04 = L_add(y22, y54);
    t07 = L_sub(y22, y54);
    t05 = L_add(y55, y23);
    t06 = L_sub(y55, y23);

    /* Post-additions */
    re[s * 3]  = L_add(t00, t04); move32();
    im[s * 3]  = L_add(t01, t05); move32();
    re[s * 11] = L_sub(t02, t06); move32();
    im[s * 11] = L_sub(t03, t07); move32();
    re[s * 19] = L_sub(t00, t04); move32();
    im[s * 19] = L_sub(t01, t05); move32();
    re[s * 27] = L_add(t02, t06); move32();
    im[s * 27] = L_add(t03, t07); move32();

    /* 5. FFT4 stage */

    /* Pre-additions */
    t00 = L_add(y08, y41);
    t02 = L_sub(y08, y41);
    t01 = L_sub(y09, y40);
    t03 = L_add(y09, y40);
    t04 = L_add(y24, y56);
    t07 = L_sub(y24, y56);
    t05 = L_add(y57, y25);
    t06 = L_sub(y57, y25);

    /* Post-additions */
    re[s * 4]  = L_add(t00, t04); move32();
    im[s * 4]  = L_add(t01, t05); move32();
    re[s * 12] = L_sub(t02, t06); move32();
    im[s * 12] = L_sub(t03, t07); move32();
    re[s * 20] = L_sub(t00, t04); move32();
    im[s * 20] = L_sub(t01, t05); move32();
    re[s * 28] = L_add(t02, t06); move32();
    im[s * 28] = L_add(t03, t07); move32();

    /* 6. FFT4 stage */

    /* Pre-additions */
    t00 = L_add(y10, y42);
    t02 = L_sub(y10, y42);
    t01 = L_add(y11, y43);
    t03 = L_sub(y11, y43);
    t04 = L_add(y26, y58);
    t07 = L_sub(y26, y58);
    t05 = L_add(y59, y27);
    t06 = L_sub(y59, y27);

    /* Post-additions */
    re[s * 5]  = L_add(t00, t04); move32();
    im[s * 5]  = L_add(t01, t05); move32();
    re[s * 13] = L_sub(t02, t06); move32();
    im[s * 13] = L_sub(t03, t07); move32();
    re[s * 21] = L_sub(t00, t04); move32();
    im[s * 21] = L_sub(t01, t05); move32();
    re[s * 29] = L_add(t02, t06); move32();
    im[s * 29] = L_add(t03, t07); move32();

    /* 7. FFT4 stage */

    /* Pre-additions */
    t00 = L_add(y12, y44);
    t02 = L_sub(y12, y44);
    t01 = L_add(y13, y45);
    t03 = L_sub(y13, y45);
    t04 = L_add(y28, y60);
    t07 = L_sub(y28, y60);
    t05 = L_add(y61, y29);
    t06 = L_sub(y61, y29);

    /* Post-additions */
    re[s * 6]  = L_add(t00, t04); move32();
    im[s * 6]  = L_add(t01, t05); move32();
    re[s * 14] = L_sub(t02, t06); move32();
    im[s * 14] = L_sub(t03, t07); move32();
    re[s * 22] = L_sub(t00, t04); move32();
    im[s * 22] = L_sub(t01, t05); move32();
    re[s * 30] = L_add(t02, t06); move32();
    im[s * 30] = L_add(t03, t07); move32();

    /* 8. FFT4 stage */

    /* Pre-additions */
    t00 = L_add(y14, y46);
    t02 = L_sub(y14, y46);
    t01 = L_add(y15, y47);
    t03 = L_sub(y15, y47);
    t04 = L_add(y30, y62);
    t07 = L_sub(y30, y62);
    t05 = L_add(y63, y31);
    t06 = L_sub(y63, y31);

    /* Post-additions */
    re[s * 7]  = L_add(t00, t04); move32();
    im[s * 7]  = L_add(t01, t05); move32();
    re[s * 15] = L_sub(t02, t06); move32();
    im[s * 15] = L_sub(t03, t07); move32();
    re[s * 23] = L_sub(t00, t04); move32();
    im[s * 23] = L_sub(t01, t05); move32();
    re[s * 31] = L_add(t02, t06); move32();
    im[s * 31] = L_add(t03, t07); move32();

    Dyn_Mem_Deluxe_Out();
}

/**
 * \brief    Function performs a complex 40-point FFT
 *           The FFT is performed inplace. The result of the FFT
 *           is scaled by SCALEFACTOR40 bits.
 *
 * \param    [i/o] re     real part
 * \param    [i/o] im     imag part
 * \param    [i  ] sx     stride real and imag part
 *
 * \return void
 */


static void fft40(Word32 *re, Word32 *im, Word16 sx, Word32 *x)
{
    Dyn_Mem_Deluxe_In(
        const Word32 *W;
        Word16        dim1, dim2;
        Counter       i, j;
        Word32        x00, x01, x02, x03, x04, x05, x06, x07, x08, x09, x10, x11, x12, x13, x14, x15;
        Word32        t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14, t15;
        Word32        s00, s01, s02, s03, s04, s05, s06, s07, s08, s09, s10, s11, s12, s13, s14, s15;
    );

    dim1 = 5; move16();
    dim2 = 8; move16();

    W = RotVector_40_32;

    FOR (i = 0; i < dim2; i++)
    {
        FOR (j = 0; j < dim1; j++)
        {
            x[2 * i * dim1 + 2 * j]     = re[sx * i + sx * j * dim2]; move32();
            x[2 * i * dim1 + 2 * j + 1] = im[sx * i + sx * j * dim2]; move32();
        }
    }

    FOR (i = 0; i < dim2; i++)
    {
        fft5(&x[i * 2 * dim1], &x[i * 2 * dim1 + 1], 2);
    }

    FOR (i = 0; i < dim1; i++)
    {
        cplxMpy4_8_1(x00, x01, x[2 * i + 2 * 0 * dim1], x[2 * i + 2 * 0 * dim1 + 1]);

        IF (i == 0)
        {
            cplxMpy4_8_1(x02, x03, x[2 * i + 2 * 1 * dim1], x[2 * i + 2 * 1 * dim1 + 1]);
            cplxMpy4_8_1(x04, x05, x[2 * i + 2 * 2 * dim1], x[2 * i + 2 * 2 * dim1 + 1]);
            cplxMpy4_8_1(x06, x07, x[2 * i + 2 * 3 * dim1], x[2 * i + 2 * 3 * dim1 + 1]);
            cplxMpy4_8_1(x08, x09, x[2 * i + 2 * 4 * dim1], x[2 * i + 2 * 4 * dim1 + 1]);
            cplxMpy4_8_1(x10, x11, x[2 * i + 2 * 5 * dim1], x[2 * i + 2 * 5 * dim1 + 1]);
            cplxMpy4_8_1(x12, x13, x[2 * i + 2 * 6 * dim1], x[2 * i + 2 * 6 * dim1 + 1]);
            cplxMpy4_8_1(x14, x15, x[2 * i + 2 * 7 * dim1], x[2 * i + 2 * 7 * dim1 + 1]);
        }
        ELSE
        {
            cplxMpy4_8_2(x02, x03, x[2 * i + 2 * 1 * dim1], x[2 * i + 2 * 1 * dim1 + 1], W[2 * (i - 1) + 0 * 2 * 4],
                         W[2 * (i - 1) + 0 * 2 * 4 + 1]);
            cplxMpy4_8_2(x04, x05, x[2 * i + 2 * 2 * dim1], x[2 * i + 2 * 2 * dim1 + 1], W[2 * (i - 1) + 1 * 2 * 4],
                         W[2 * (i - 1) + 1 * 2 * 4 + 1]);
            cplxMpy4_8_2(x06, x07, x[2 * i + 2 * 3 * dim1], x[2 * i + 2 * 3 * dim1 + 1], W[2 * (i - 1) + 2 * 2 * 4],
                         W[2 * (i - 1) + 2 * 2 * 4 + 1]);
            cplxMpy4_8_2(x08, x09, x[2 * i + 2 * 4 * dim1], x[2 * i + 2 * 4 * dim1 + 1], W[2 * (i - 1) + 3 * 2 * 4],
                         W[2 * (i - 1) + 3 * 2 * 4 + 1]);
            cplxMpy4_8_2(x10, x11, x[2 * i + 2 * 5 * dim1], x[2 * i + 2 * 5 * dim1 + 1], W[2 * (i - 1) + 4 * 2 * 4],
                         W[2 * (i - 1) + 4 * 2 * 4 + 1]);
            cplxMpy4_8_2(x12, x13, x[2 * i + 2 * 6 * dim1], x[2 * i + 2 * 6 * dim1 + 1], W[2 * (i - 1) + 5 * 2 * 4],
                         W[2 * (i - 1) + 5 * 2 * 4 + 1]);
            cplxMpy4_8_2(x14, x15, x[2 * i + 2 * 7 * dim1], x[2 * i + 2 * 7 * dim1 + 1], W[2 * (i - 1) + 6 * 2 * 4],
                         W[2 * (i - 1) + 6 * 2 * 4 + 1]);
        }

        t00 = L_shr_pos(L_add(x00, x08), SCALEFACTORN2 - 1);
        t02 = L_shr_pos(L_sub(x00, x08), SCALEFACTORN2 - 1);
        t01 = L_shr_pos(L_add(x01, x09), SCALEFACTORN2 - 1);
        t03 = L_shr_pos(L_sub(x01, x09), SCALEFACTORN2 - 1);
        t04 = L_shr_pos(L_add(x02, x10), SCALEFACTORN2 - 1);
        t06 = L_sub(x02, x10);
        t05 = L_shr_pos(L_add(x03, x11), SCALEFACTORN2 - 1);
        t07 = L_sub(x03, x11);
        t08 = L_shr_pos(L_add(x04, x12), SCALEFACTORN2 - 1);
        t10 = L_shr_pos(L_sub(x04, x12), SCALEFACTORN2 - 1);
        t09 = L_shr_pos(L_add(x05, x13), SCALEFACTORN2 - 1);
        t11 = L_shr_pos(L_sub(x05, x13), SCALEFACTORN2 - 1);
        t12 = L_shr_pos(L_add(x06, x14), SCALEFACTORN2 - 1);
        t14 = L_sub(x06, x14);
        t13 = L_shr_pos(L_add(x07, x15), SCALEFACTORN2 - 1);
        t15 = L_sub(x07, x15);

        s00 = L_add(t00, t08);
        s04 = L_sub(t00, t08);
        s01 = L_add(t01, t09);
        s05 = L_sub(t01, t09);
        s08 = L_sub(t02, t11);
        s10 = L_add(t02, t11);
        s09 = L_add(t03, t10);
        s11 = L_sub(t03, t10);
        s02 = L_add(t04, t12);
        s07 = L_sub(t04, t12);
        s03 = L_add(t05, t13);
        s06 = L_sub(t13, t05);

        t01 = L_shr_pos(L_add(t06, t14), SCALEFACTORN2 - 1);
        t02 = L_shr_pos(L_sub(t06, t14), SCALEFACTORN2 - 1);
        t00 = L_shr_pos(L_add(t07, t15), SCALEFACTORN2 - 1);
        t03 = L_shr_pos(L_sub(t07, t15), SCALEFACTORN2 - 1);

        s12 = Mpy_32_32(L_add(t00, t02), C81_32);
        s14 = Mpy_32_32(L_sub(t00, t02), C81_32);
        s13 = Mpy_32_32(L_sub(t03, t01), C81_32);
        s15 = Mpy_32_32(L_add(t01, t03), C82_32);

        re[sx * i + sx * 0 * dim1] = L_add(s00, s02); move32();
        im[sx * i + sx * 0 * dim1] = L_add(s01, s03); move32();
        re[sx * i + sx * 1 * dim1] = L_add(s10, s12); move32();
        im[sx * i + sx * 1 * dim1] = L_add(s11, s13); move32();
        re[sx * i + sx * 2 * dim1] = L_sub(s04, s06); move32();
        im[sx * i + sx * 2 * dim1] = L_sub(s05, s07); move32();
        re[sx * i + sx * 3 * dim1] = L_add(s08, s14); move32();
        im[sx * i + sx * 3 * dim1] = L_add(s09, s15); move32();
        re[sx * i + sx * 4 * dim1] = L_sub(s00, s02); move32();
        im[sx * i + sx * 4 * dim1] = L_sub(s01, s03); move32();
        re[sx * i + sx * 5 * dim1] = L_sub(s10, s12); move32();
        im[sx * i + sx * 5 * dim1] = L_sub(s11, s13); move32();
        re[sx * i + sx * 6 * dim1] = L_add(s04, s06); move32();
        im[sx * i + sx * 6 * dim1] = L_add(s05, s07); move32();
        re[sx * i + sx * 7 * dim1] = L_sub(s08, s14); move32();
        im[sx * i + sx * 7 * dim1] = L_sub(s09, s15); move32();
    }

    Dyn_Mem_Deluxe_Out();
}

/**
 * \brief Combined FFT
 *
 * \param    [i/o] re     real part
 * \param    [i/o] im     imag part
 * \param    [i  ] W      rotation factor
 * \param    [i  ] dim1   length of fft1
 * \param    [i  ] dim2   length of fft2
 * \param    [i  ] sx     stride real and imag part
 * \param    [i  ] sc     stride phase rotation coefficients
 * \param    [tmp] x      32-bit workbuffer of length=2*len
 * \param    [i  ] Woff   offset for addressing the rotation vector table
 *
 * \return void
 */


static void fftN2(Word32 *re, Word32 *im, const Word16 *W, Word16 dim1, Word16 dim2, Word16 sx, Word16 sc, Word32 *x,
                  Word16 Woff)
{
    Dyn_Mem_Deluxe_In(
        Counter i, j;
    );


    FOR (i = 0; i < dim2; i++)
    {
        FOR (j = 0; j < dim1; j++)
        {
            x[2 * i * dim1 + 2 * j]     = re[sx * i + sx * j * dim2]; move32();
            x[2 * i * dim1 + 2 * j + 1] = im[sx * i + sx * j * dim2]; move32();
        }
    }

    SWITCH (dim1)
    {
    case 4:
        FOR (i = 0; i < dim2; i++)
        {
            fft4(&x[i * 2 * dim1]);
        }
        BREAK;
    case 8:
        FOR (i = 0; i < dim2; i++)
        {
            fft8(&x[i * 2 * dim1], &x[i * 2 * dim1 + 1], 2);
        }
        BREAK;
    case 10:
        FOR (i = 0; i < dim2; i++)
        {
            fft10(&x[i * 2 * dim1], &x[i * 2 * dim1 + 1], 2);
        }
        BREAK;
    case 15:
        FOR (i = 0; i < dim2; i++)
        {
            fft15(&x[i * 2 * dim1], &x[i * 2 * dim1 + 1], 2);
        }
        BREAK;
    case 16:
        FOR (i = 0; i < dim2; i++)
        {
            fft16(&x[i * 2 * dim1], &x[i * 2 * dim1 + 1], 2);
        }
        BREAK;
    case 20:
        FOR (i = 0; i < dim2; i++)
        {
            fft20(&x[i * 2 * dim1], &x[i * 2 * dim1 + 1], 2);
        }
        BREAK;
    case 30:
        FOR (i = 0; i < dim2; i++)
        {
            fft30(&x[i * 2 * dim1], &x[i * 2 * dim1 + 1], 2);
        }
        BREAK;
    case 32:
        FOR (i = 0; i < dim2; i++)
        {
            fft32(&x[i * 2 * dim1], &x[i * 2 * dim1 + 1], 2);
        }
        BREAK;
    default: ASSERT(0);
    }

    SWITCH (dim2)
    {
    case 4:
    {
        Word32 x00, x01, x02, x03, x04, x05, x06, x07;
        Word32 t00, t01, t02, t03, t04, t05, t06, t07;

        j = add(8, 0);
        FOR (i = 0; i < dim1; i++)
        {
            cplxMpy4_4_1(x00, x01, x[2 * i + 2 * 0 * dim1], x[2 * i + 2 * 0 * dim1 + 1]);
            IF (i == 0)
            {
                cplxMpy4_4_1(x02, x03, x[2 * i + 2 * 1 * dim1], x[2 * i + 2 * 1 * dim1 + 1]);
                cplxMpy4_4_1(x04, x05, x[2 * i + 2 * 2 * dim1], x[2 * i + 2 * 2 * dim1 + 1]);
                cplxMpy4_4_1(x06, x07, x[2 * i + 2 * 3 * dim1], x[2 * i + 2 * 3 * dim1 + 1]);
            }
            ELSE
            {
                cplxMpy4_4_0(x02, x03, x[2 * i + 2 * 1 * dim1], x[2 * i + 2 * 1 * dim1 + 1],
                             W[sc * i + j * 1 * dim1 - Woff], W[sc * i + j * 1 * dim1 + 1 - Woff]);
                cplxMpy4_4_0(x04, x05, x[2 * i + 2 * 2 * dim1], x[2 * i + 2 * 2 * dim1 + 1],
                             W[sc * i + j * 2 * dim1 - Woff], W[sc * i + j * 2 * dim1 + 1 - Woff]);
                cplxMpy4_4_0(x06, x07, x[2 * i + 2 * 3 * dim1], x[2 * i + 2 * 3 * dim1 + 1],
                             W[sc * i + j * 3 * dim1 - Woff], W[sc * i + j * 3 * dim1 + 1 - Woff]);
            }

            t00 = L_add(x00, x04);
            t02 = L_sub(x00, x04);
            t01 = L_add(x01, x05);
            t03 = L_sub(x01, x05);
            t04 = L_add(x02, x06);
            t07 = L_sub(x02, x06);
            t05 = L_add(x07, x03);
            t06 = L_sub(x07, x03);

            re[sx * i + sx * 0 * dim1] = L_add(t00, t04); move32();
            im[sx * i + sx * 0 * dim1] = L_add(t01, t05); move32();
            re[sx * i + sx * 1 * dim1] = L_sub(t02, t06); move32();
            im[sx * i + sx * 1 * dim1] = L_sub(t03, t07); move32();
            re[sx * i + sx * 2 * dim1] = L_sub(t00, t04); move32();
            im[sx * i + sx * 2 * dim1] = L_sub(t01, t05); move32();
            re[sx * i + sx * 3 * dim1] = L_add(t02, t06); move32();
            im[sx * i + sx * 3 * dim1] = L_add(t03, t07); move32();
        }
        BREAK;
    }
    case 8:
    {
        Word32 x00, x01, x02, x03, x04, x05, x06, x07, x08, x09, x10, x11, x12, x13, x14, x15;
        Word32 t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14, t15;
        Word32 s00, s01, s02, s03, s04, s05, s06, s07, s08, s09, s10, s11, s12, s13, s14, s15;

        FOR (i = 0; i < dim1; i++)
        {
            cplxMpy4_8_1(x00, x01, x[2 * i + 2 * 0 * dim1], x[2 * i + 2 * 0 * dim1 + 1]);
            IF (i == 0)
            {
                cplxMpy4_8_1(x02, x03, x[2 * i + 2 * 1 * dim1], x[2 * i + 2 * 1 * dim1 + 1]);
                cplxMpy4_8_1(x04, x05, x[2 * i + 2 * 2 * dim1], x[2 * i + 2 * 2 * dim1 + 1]);
                cplxMpy4_8_1(x06, x07, x[2 * i + 2 * 3 * dim1], x[2 * i + 2 * 3 * dim1 + 1]);
                cplxMpy4_8_1(x08, x09, x[2 * i + 2 * 4 * dim1], x[2 * i + 2 * 4 * dim1 + 1]);
                cplxMpy4_8_1(x10, x11, x[2 * i + 2 * 5 * dim1], x[2 * i + 2 * 5 * dim1 + 1]);
                cplxMpy4_8_1(x12, x13, x[2 * i + 2 * 6 * dim1], x[2 * i + 2 * 6 * dim1 + 1]);
                cplxMpy4_8_1(x14, x15, x[2 * i + 2 * 7 * dim1], x[2 * i + 2 * 7 * dim1 + 1]);
            }
            ELSE
            {
                cplxMpy4_8_0(x02, x03, x[2 * i + 2 * 1 * dim1], x[2 * i + 2 * 1 * dim1 + 1],
                             W[sc * i + sc * 1 * dim1 - Woff], W[sc * i + sc * 1 * dim1 + 1 - Woff]);
                cplxMpy4_8_0(x04, x05, x[2 * i + 2 * 2 * dim1], x[2 * i + 2 * 2 * dim1 + 1],
                             W[sc * i + sc * 2 * dim1 - Woff], W[sc * i + sc * 2 * dim1 + 1 - Woff]);
                cplxMpy4_8_0(x06, x07, x[2 * i + 2 * 3 * dim1], x[2 * i + 2 * 3 * dim1 + 1],
                             W[sc * i + sc * 3 * dim1 - Woff], W[sc * i + sc * 3 * dim1 + 1 - Woff]);
                cplxMpy4_8_0(x08, x09, x[2 * i + 2 * 4 * dim1], x[2 * i + 2 * 4 * dim1 + 1],
                             W[sc * i + sc * 4 * dim1 - Woff], W[sc * i + sc * 4 * dim1 + 1 - Woff]);
                cplxMpy4_8_0(x10, x11, x[2 * i + 2 * 5 * dim1], x[2 * i + 2 * 5 * dim1 + 1],
                             W[sc * i + sc * 5 * dim1 - Woff], W[sc * i + sc * 5 * dim1 + 1 - Woff]);
                cplxMpy4_8_0(x12, x13, x[2 * i + 2 * 6 * dim1], x[2 * i + 2 * 6 * dim1 + 1],
                             W[sc * i + sc * 6 * dim1 - Woff], W[sc * i + sc * 6 * dim1 + 1 - Woff]);
                cplxMpy4_8_0(x14, x15, x[2 * i + 2 * 7 * dim1], x[2 * i + 2 * 7 * dim1 + 1],
                             W[sc * i + sc * 7 * dim1 - Woff], W[sc * i + sc * 7 * dim1 + 1 - Woff]);
            }

            t00 = L_shr_pos(L_add(x00, x08), SCALEFACTORN2 - 1);
            t02 = L_shr_pos(L_sub(x00, x08), SCALEFACTORN2 - 1);
            t01 = L_shr_pos(L_add(x01, x09), SCALEFACTORN2 - 1);
            t03 = L_shr_pos(L_sub(x01, x09), SCALEFACTORN2 - 1);
            t04 = L_shr_pos(L_add(x02, x10), SCALEFACTORN2 - 1);
            t06 = L_sub(x02, x10);
            t05 = L_shr_pos(L_add(x03, x11), SCALEFACTORN2 - 1);
            t07 = L_sub(x03, x11);
            t08 = L_shr_pos(L_add(x04, x12), SCALEFACTORN2 - 1);
            t10 = L_shr_pos(L_sub(x04, x12), SCALEFACTORN2 - 1);
            t09 = L_shr_pos(L_add(x05, x13), SCALEFACTORN2 - 1);
            t11 = L_shr_pos(L_sub(x05, x13), SCALEFACTORN2 - 1);
            t12 = L_shr_pos(L_add(x06, x14), SCALEFACTORN2 - 1);
            t14 = L_sub(x06, x14);
            t13 = L_shr_pos(L_add(x07, x15), SCALEFACTORN2 - 1);
            t15 = L_sub(x07, x15);

            s00 = L_add(t00, t08);
            s04 = L_sub(t00, t08);
            s01 = L_add(t01, t09);
            s05 = L_sub(t01, t09);
            s08 = L_sub(t02, t11);
            s10 = L_add(t02, t11);
            s09 = L_add(t03, t10);
            s11 = L_sub(t03, t10);
            s02 = L_add(t04, t12);
            s07 = L_sub(t04, t12);
            s03 = L_add(t05, t13);
            s06 = L_sub(t13, t05);

            t01 = L_shr_pos(L_add(t06, t14), SCALEFACTORN2 - 1);
            t02 = L_shr_pos(L_sub(t06, t14), SCALEFACTORN2 - 1);
            t00 = L_shr_pos(L_add(t07, t15), SCALEFACTORN2 - 1);
            t03 = L_shr_pos(L_sub(t07, t15), SCALEFACTORN2 - 1);

            s12 = Mpy_32_xx(L_add(t00, t02), C81);
            s14 = Mpy_32_xx(L_sub(t00, t02), C81);
            s13 = Mpy_32_xx(L_sub(t03, t01), C81);
            s15 = Mpy_32_xx(L_add(t01, t03), C82);

            re[sx * i + sx * 0 * dim1] = L_add(s00, s02); move32();
            im[sx * i + sx * 0 * dim1] = L_add(s01, s03); move32();
            re[sx * i + sx * 1 * dim1] = L_add(s10, s12); move32();
            im[sx * i + sx * 1 * dim1] = L_add(s11, s13); move32();
            re[sx * i + sx * 2 * dim1] = L_sub(s04, s06); move32();
            im[sx * i + sx * 2 * dim1] = L_sub(s05, s07); move32();
            re[sx * i + sx * 3 * dim1] = L_add(s08, s14); move32();
            im[sx * i + sx * 3 * dim1] = L_add(s09, s15); move32();
            re[sx * i + sx * 4 * dim1] = L_sub(s00, s02); move32();
            im[sx * i + sx * 4 * dim1] = L_sub(s01, s03); move32();
            re[sx * i + sx * 5 * dim1] = L_sub(s10, s12); move32();
            im[sx * i + sx * 5 * dim1] = L_sub(s11, s13); move32();
            re[sx * i + sx * 6 * dim1] = L_add(s04, s06); move32();
            im[sx * i + sx * 6 * dim1] = L_add(s05, s07); move32();
            re[sx * i + sx * 7 * dim1] = L_sub(s08, s14); move32();
            im[sx * i + sx * 7 * dim1] = L_sub(s09, s15); move32();
        }
        BREAK;
    }

    case 12:
    {
        Word32 y[2 * 20];
        FOR (j = 0; j < dim2; j++)
        {
            cplxMpy4_12_1(y[2 * j], y[2 * j + 1], x[2 * 0 + 2 * j * dim1], x[2 * 0 + 2 * j * dim1 + 1]);
        }
        fft12(y);
        FOR (j = 0; j < dim2; j++)
        {
            re[sx * 0 + sx * j * dim1] = y[2 * j];     move32();
            im[sx * 0 + sx * j * dim1] = y[2 * j + 1]; move32();
        }

        FOR (i = 1; i < dim1; i++)
        {
            cplxMpy4_12_1(y[2 * (0 + 0)], y[2 * (0 + 0) + 1], x[2 * i + 2 * (0 + 0) * dim1],
                          x[2 * i + 2 * (0 + 0) * dim1 + 1]);
            cplxMpy4_12_0(y[2 * (0 + 1)], y[2 * (0 + 1) + 1], x[2 * i + 2 * (0 + 1) * dim1],
                          x[2 * i + 2 * (0 + 1) * dim1 + 1], W[sc * i + sc * (0 + 1) * dim1 - Woff],
                          W[sc * i + sc * (0 + 1) * dim1 + 1 - Woff]);
            FOR (j = 2; j < dim2; j = j + 2)
            {
                cplxMpy4_12_0(y[2 * (j + 0)], y[2 * (j + 0) + 1], x[2 * i + 2 * (j + 0) * dim1],
                              x[2 * i + 2 * (j + 0) * dim1 + 1], W[sc * i + sc * (j + 0) * dim1 - Woff],
                              W[sc * i + sc * (j + 0) * dim1 + 1 - Woff]);
                cplxMpy4_12_0(y[2 * (j + 1)], y[2 * (j + 1) + 1], x[2 * i + 2 * (j + 1) * dim1],
                              x[2 * i + 2 * (j + 1) * dim1 + 1], W[sc * i + sc * (j + 1) * dim1 - Woff],
                              W[sc * i + sc * (j + 1) * dim1 + 1 - Woff]);
            }
            fft12(y);
            FOR (j = 0; j < dim2; j++)
            {
                re[sx * i + sx * j * dim1] = y[2 * j];     move32();
                im[sx * i + sx * j * dim1] = y[2 * j + 1]; move32();
            }
        }
        BREAK;
    }
    default: ASSERT(0);
    }

    Dyn_Mem_Deluxe_Out();
}

/**
 * \brief Complex valued FFT
 *
 * \param    [i/o] re          real part
 * \param    [i/o] im          imag part
 * \param    [i  ] sizeOfFft   length of fft
 * \param    [i  ] s           stride real and imag part
 * \param    [i  ] scale       scalefactor
 *
 * \return void
 */


/* x is the scratch buffer */
void BASOP_cfft(Word32 *re, Word32 *im, Word16 length, Word16 s, Word16 *scale, Word32 *x)
{
    SWITCH (length)
    {
    case 10:
        fft10(re, im, s);
        *scale = add(*scale, SCALEFACTOR10); move16();
        BREAK;
    case 16:
        fft16(re, im, s);
        *scale = add(*scale, SCALEFACTOR16); move16();
        BREAK;
    case 20:
        fft20(re, im, s);
        *scale = add(*scale, SCALEFACTOR20); move16();
        BREAK;
    case 30:
        fft30(re, im, s);
        *scale = add(*scale, SCALEFACTOR30); move16();
        BREAK;
    case 32:
        fft32(re, im, s);
        *scale = add(*scale, SCALEFACTOR32); move16();
        BREAK;
    case 40:
        fft40(re, im, s, x);
        *scale = add(*scale, SCALEFACTOR40); move16();
        BREAK;
    case 48:
        fftN2(re, im, RotVector_32_12, 4, 12, s, 16, x, 64);
        *scale = add(*scale, SCALEFACTOR48); move16();
        BREAK;
    case 60:
        fftN2(re, im, RotVector_480, 15, 4, s, 4, x, 60);
        *scale = add(*scale, SCALEFACTOR60); move16();
        BREAK;
    case 64:
        fftN2(re, im, RotVector_32_8, 8, 8, s, 8, x, 64);
        *scale = add(*scale, SCALEFACTOR64); move16();
        BREAK;
    case 80:
        fftN2(re, im, RotVector_320, 10, 8, s, 4, x, 40);
        *scale = add(*scale, SCALEFACTOR80); move16();
        BREAK;
    case 96:
        fftN2(re, im, RotVector_32_12, 8, 12, s, 8, x, 64);
        *scale = add(*scale, SCALEFACTOR96); move16();
        BREAK;
    case 120:
        fftN2(re, im, RotVector_480, 15, 8, s, 4, x, 60);
        *scale = add(*scale, SCALEFACTOR120); move16();
        BREAK;
    case 128:
        fftN2(re, im, RotVector_32_8, 16, 8, s, 4, x, 64);
        *scale = add(*scale, SCALEFACTOR128); move16();
        BREAK;
    case 160:
        fftN2(re, im, RotVector_320, 20, 8, s, 2, x, 40);
        *scale = add(*scale, SCALEFACTOR160); move16();
        BREAK;
    case 192:
        fftN2(re, im, RotVector_32_12, 16, 12, s, 4, x, 64);
        *scale = add(*scale, SCALEFACTOR192); move16();
        BREAK;
    case 240:
        fftN2(re, im, RotVector_480, 30, 8, s, 2, x, 60);
        *scale = add(*scale, SCALEFACTOR240); move16();
        BREAK;
    case 256:
        fftN2(re, im, RotVector_32_8, 32, 8, s, 2, x, 64);
        *scale = add(*scale, SCALEFACTOR256); move16();
        BREAK;
    case 384:
        fftN2(re, im, RotVector_32_12, 32, 12, s, 2, x, 64);
        *scale = add(*scale, SCALEFACTOR384); move16();
        BREAK;
    default: ASSERT(0);
    }
}

#endif /* !USE_KISS_FFT */

#define RFFT_TWIDDLE1(x, t1, t2, t3, t4, w1, w2, xb0, xb1, xt0, xt1)                                                   \
    do                                                                                                                 \
    {                                                                                                                  \
        xb0 = L_shr_pos(x[2 * i + 0], 2);                                                                              \
        xb1 = L_shr_pos(x[2 * i + 1], 2);                                                                              \
        xt0 = L_shr_pos(x[sizeOfFft - 2 * i + 0], 2);                                                                  \
        xt1 = L_shr_pos(x[sizeOfFft - 2 * i + 1], 2);                                                                  \
        t1  = L_sub(xb0, xt0);                                                                                         \
        t2  = L_add(xb1, xt1);                                                                                         \
        t3  = L_sub(Mpy_32_32(t1, w1), Mpy_32_32(t2, w2));                                                             \
        t4  = L_add(Mpy_32_32(t1, w2), Mpy_32_32(t2, w1));                                                             \
        t1  = L_add(xb0, xt0);                                                                                         \
        t2  = L_sub(xb1, xt1);                                                                                         \
    } while (0)

#define RFFT_TWIDDLE2(x, t1, t2, t3, t4, w1, w2, xb0, xb1, xt0, xt1)                                                   \
    do                                                                                                                 \
    {                                                                                                                  \
        xb0 = L_shr_pos(x[2 * i + 0], 2);                                                                              \
        xb1 = L_shr_pos(x[2 * i + 1], 2);                                                                              \
        xt0 = L_shr_pos(x[sizeOfFft - 2 * i + 0], 2);                                                                  \
        xt1 = L_shr_pos(x[sizeOfFft - 2 * i + 1], 2);                                                                  \
        t1  = L_sub(xb0, xt0);                                                                                         \
        t2  = L_add(xb1, xt1);                                                                                         \
        t3  = L_add(Mpy_32_32(t1, w1), Mpy_32_32(t2, w2));                                                             \
        t4  = L_sub(Mpy_32_32(t2, w1), Mpy_32_32(t1, w2));                                                             \
        t1  = L_add(xb0, xt0);                                                                                         \
        t2  = L_sub(xb1, xt1);                                                                                         \
    } while (0)



static const Word32 *rfft_twid(int size)
{
    SWITCH (size)
    {
    case 32: return RealFFT32_twid;
    case 40: return RealFFT40_twid;
    case 64: return RealFFT64_twid;
    case 80: return RealFFT80_twid;
    case 96: return RealFFT96_twid;
    case 128: return RealFFT128_twid;
    case 192: return RealFFT192_twid;
    case 256: return RealFFT256_twid;
    case 384: return RealFFT384_twid;
    case 512: return RealFFT512_twid;
    case 768: return RealFFT768_twid;
    default: ASSERT(0);
    }
    return NULL;
}


void BASOP_rfftN(Word32 *x, Word16 sizeOfFft, Word16 *scale, Word8 *scratchBuffer)
{
    Dyn_Mem_Deluxe_In(
        Counter       i;
        Word16        sizeOfFft2, sizeOfFft4, sizeOfFft8;
        Word32        t1, t2, t3, t4, xb0, xb1, xt0, xt1;
        Word32 *      workBuffer;
        const Word32 *w32;
    );

    workBuffer = (Word32 *)scratchAlign(scratchBuffer, 0); /* Size = 4 * sizeOfFft */
    w32        = rfft_twid(sizeOfFft);

    sizeOfFft2 = shr_pos(sizeOfFft, 1);
    sizeOfFft4 = shr_pos(sizeOfFft, 2);
    sizeOfFft8 = shr_pos(sizeOfFft, 3);

    BASOP_cfft(&x[0], &x[1], sizeOfFft2, 2, scale, workBuffer);

    xb0  = L_shr_pos(x[0], 1);
    xb1  = L_shr_pos(x[1], 1);
    x[0] = L_add(xb0, xb1); move32();
    x[1] = L_sub(xb0, xb1); move32();

    FOR (i = 1; i < sizeOfFft8; i++)
    {
        RFFT_TWIDDLE1(x, t1, t2, t3, t4, w32[2 * i + 1], w32[2 * i], xb0, xb1, xt0, xt1);
        x[2 * i]                 = L_sub(t1, t3);           move32();
        x[2 * i + 1]             = L_sub(t2, t4);           move32();
        x[sizeOfFft - 2 * i]     = L_add(t1, t3);           move32();
        x[sizeOfFft - 2 * i + 1] = L_negate(L_add(t2, t4)); move32();
    }

    FOR (i = sizeOfFft8; i < sizeOfFft4; i++)
    {
        RFFT_TWIDDLE1(x, t1, t2, t3, t4, w32[(2 * sizeOfFft4 - 2 * i)], w32[(2 * sizeOfFft4 - 2 * i + 1)], xb0, xb1,
                      xt0, xt1);
        x[2 * i]                 = L_sub(t1, t3);           move32();
        x[2 * i + 1]             = L_sub(t2, t4);           move32();
        x[sizeOfFft - 2 * i]     = L_add(t1, t3);           move32();
        x[sizeOfFft - 2 * i + 1] = L_negate(L_add(t2, t4)); move32();
    }

    x[sizeOfFft - 2 * i]     = L_shr_pos(x[2 * i + 0], 1);           move32();
    x[sizeOfFft - 2 * i + 1] = L_negate(L_shr_pos(x[2 * i + 1], 1)); move32();

    *scale = add(*scale, 1); move16();

    Dyn_Mem_Deluxe_Out();
}



void BASOP_irfftN(Word32 *x, Word16 sizeOfFft, Word16 *scale, Word8 *scratchBuffer)
{
    Dyn_Mem_Deluxe_In(
        Word16        sizeOfFft2, sizeOfFft4, sizeOfFft8;
        Word32        t1, t2, t3, t4, xb0, xb1, xt0, xt1;
        Word32 *      workBuffer;
        const Word32 *w32;
        Counter       i;
    );

    workBuffer = (Word32 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * BASOP_CFFT_MAX_LENGTH */

    w32 = rfft_twid(sizeOfFft);

    sizeOfFft2 = shr_pos(sizeOfFft, 1);
    sizeOfFft4 = shr_pos(sizeOfFft, 2);
    sizeOfFft8 = shr_pos(sizeOfFft, 3);

    xb0  = L_shr_pos(x[0], 2);
    xb1  = L_shr_pos(x[1], 2);
    x[0] = L_add(xb0, xb1); move32();
    x[1] = L_sub(xb1, xb0); move32();

    FOR (i = 1; i < sizeOfFft8; i++)
    {
        RFFT_TWIDDLE2(x, t1, t2, t3, t4, w32[2 * i + 1], w32[2 * i], xb0, xb1, xt0, xt1);
        x[2 * i]                 = L_sub(t1, t3); move32();
        x[2 * i + 1]             = L_sub(t4, t2); move32();
        x[sizeOfFft - 2 * i]     = L_add(t1, t3); move32();
        x[sizeOfFft - 2 * i + 1] = L_add(t2, t4); move32();
    }

    FOR (i = sizeOfFft8; i < sizeOfFft4; i++)
    {
        RFFT_TWIDDLE2(x, t1, t2, t3, t4, w32[(2 * sizeOfFft4 - 2 * i)], w32[(2 * sizeOfFft4 - 2 * i + 1)], xb0, xb1,
                      xt0, xt1);
        x[2 * i]                 = L_sub(t1, t3); move32();
        x[2 * i + 1]             = L_sub(t4, t2); move32();
        x[sizeOfFft - 2 * i]     = L_add(t1, t3); move32();
        x[sizeOfFft - 2 * i + 1] = L_add(t2, t4); move32();
    }

    x[sizeOfFft - 2 * i]     = L_shr_pos(x[2 * i + 0], 1); move32();
    x[sizeOfFft - 2 * i + 1] = L_shr_pos(x[2 * i + 1], 1); move32();

    BASOP_cfft(&x[0], &x[1], sizeOfFft2, 2, scale, workBuffer);

    /* If you want BASOP_irfft to be inverse to BASOP_rfft then the result needs
     * to be normalised by sizeOfFft */
    FOR (i = 0; i < sizeOfFft2; i++)
    {
        x[2 * i + 1] = L_negate(x[2 * i + 1]); move32();
    }

    *scale = add(*scale, 2); move16();

    Dyn_Mem_Deluxe_Out();
}


#ifdef USE_KISS_FFT
/* experimental fft backend using kiss fft
 * download from https://sourceforge.net/projects/kissfft/
 * extract kiss_fft130.zip in directory which contains fft.c */
#define FIXED_POINT 32
#include "kiss_fft130/kiss_fft.c"

typedef long long Word64;
/* cache fft configs for reuse. kfft_fft_cfg kfft_scalefactors and kfft_idx must be changed in conjunction */
static kiss_fft_cfg kfft_fft_cfg[22];

static const int kfft_scalefactors[22] = {
    SCALEFACTOR10,
    SCALEFACTOR16,
    SCALEFACTOR20,
    SCALEFACTOR30,
    SCALEFACTOR32,
    SCALEFACTOR40,
    SCALEFACTOR48,
    SCALEFACTOR60,
    SCALEFACTOR64,
    SCALEFACTOR80,
    SCALEFACTOR90,
    SCALEFACTOR96,
    SCALEFACTOR120,
    SCALEFACTOR128,
    SCALEFACTOR160,
    SCALEFACTOR180,
    SCALEFACTOR192,
    SCALEFACTOR240,
    SCALEFACTOR256,
    SCALEFACTOR384,
    -1,
    -1,
};


static int kfft_idx(int len)
{
    switch (len)
    {
    case 10: return 0;
    case 16: return 1;
    case 20: return 2;
    case 30: return 3;
    case 32: return 4;
    case 40: return 5;
    case 48: return 6;
    case 60: return 7;
    case 64: return 8;
    case 80: return 9;
    case 90: return 10;
    case 96: return 11;
    case 120: return 12;
    case 128: return 13;
    case 160: return 14;
    case 180: return 15;
    case 192: return 16;
    case 240: return 17;
    case 256: return 18;
    case 384: return 19;
    case 512: return 20;
    case 768: return 21;
    default: ASSERT(0);
    }
    return 0;
}

/* free all kiss fft handles. */
static void kfft_free(void)
{
    const int numcfg = sizeof(kfft_fft_cfg) / sizeof(*kfft_fft_cfg);
    for (int idx = 0; idx < numcfg; idx++)
    {
        kiss_fft_free(kfft_fft_cfg[idx]);
        kfft_fft_cfg[idx] = NULL;
    }
}


static void kfft_fft(Word32 *restrict re, Word32 *restrict im, int len, int stride, int rshift,
                     Word32 *restrict scratch)
{
    int idx = kfft_idx(len);
    /* create config for length if it doesn't exist, not thread safe! */
    if (!kfft_fft_cfg[idx])
    {
        static int first_run = 1;

        kfft_fft_cfg[idx] = kiss_fft_alloc(len, 0, NULL, NULL);
        if (first_run)
        {
            atexit(kfft_free);
            first_run = 0;
        }
    }
    /* interleave input if required, else do nothing */
    kiss_fft_cpx *buf = NULL;
    if (re + 1 == im && stride == 2)
    {
        buf = (kiss_fft_cpx *)re;
    }
    else
    {
        buf = (kiss_fft_cpx *)scratch;
        for (int i = 0; i < len; i++)
        {
            buf[i].r = re[i * stride];
            buf[i].i = im[i * stride];
        }
    }

    kiss_fft(kfft_fft_cfg[idx], buf, buf);

    /* undo normalisation for compability */
    for (int i = 0; i < len; i++)
    {
        re[i * stride] = (Word32)(((Word64)buf[i].r * len) >> rshift);
        im[i * stride] = (Word32)(((Word64)buf[i].i * len) >> rshift);
    }
}


void BASOP_cfft(Word32 *re, Word32 *im, Word16 sizeOfFft, Word16 s, Word16 *scale, Word32 *x)
{
    int rshift = kfft_scalefactors[kfft_idx(sizeOfFft)];
    ASSERT(rshift >= 0);
    kfft_fft(re, im, sizeOfFft, s, rshift, x);
    *scale += rshift;
}


void fft16(Word32 *re, Word32 *im, Word16 s)
{
    Word32 scratch[32];
    kfft_fft(re, im, 16, s, 0, scratch);
}


#endif /* USE_KISS_FFT */
