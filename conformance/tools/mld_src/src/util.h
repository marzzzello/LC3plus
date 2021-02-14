/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef UTIL_H
#define UTIL_H

#include "shims.h"

// number of elements in array
#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

// can't use complex.h, MSVC is still stuck in 1989 :(
typedef struct {
    double r;
    double i;
} Complex;

// complex constructor
static inline Complex cmplx(double r, double i) { return (Complex){r, i}; }

// complex a + b
static inline Complex cadd(Complex a, Complex b) { return cmplx(a.r + b.r, a.i + b.i); }

// complex a * b
static inline Complex cmul(Complex a, Complex b) { return cmplx(a.r * b.r - a.i * b.i, a.i * b.r + a.r * b.i); }

// complex eᶦˣ
static inline Complex cexpi(double x) { return cmplx(cos(x), sin(x)); }

// complex -x
static inline Complex cneg(Complex x) { return cmplx(-x.r, -x.i); }

// min without side effects
static inline int imin(int a, int b) { return a < b ? a : b; }

// 10ˣ
static inline double exp10(double x) { return pow(10, x); }

// x²
static inline double sqr(double x) { return x * x; }

// atof with error checking
static inline bool atof_c(const char* str, double* value)
{
    char* end = NULL;
    *value    = str ? strtod(str, &end) : 0;
    return str && *end == 0;
}

#endif
