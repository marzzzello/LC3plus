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

#include "clib.h"
#include "math.h"

#ifdef _MSC_VER
/* strcasecmp is not available on visual studio */
static LC3_INT strcasecmp(const char* a, const char* b) {
  return _stricmp(a,b);
}
#endif

/* restrict is not available on visual studio */
#ifdef _MSC_VER
#define restrict __restrict
#endif

/* number of elements in array */
#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

/* min max with no side effects */
static inline LC3_INT imin(LC3_INT a, LC3_INT b) { return a < b ? a : b; }
static inline LC3_INT imax(LC3_INT a, LC3_INT b) { return a > b ? a : b; }

/* restrict x to range [min, max] */
static inline LC3_INT iclamp(LC3_INT min, LC3_INT x, LC3_INT max) {
  return x < min ? min : x > max ? max : x;
}
static inline double fcmamp(double min, double x, double max) {
  return x < min ? min : x > max ? max : x;
}
static inline LC3_FLOAT fclampf(LC3_FLOAT min, LC3_FLOAT x, LC3_FLOAT max) {
  return x < min ? min : x > max ? max : x;
}

/* x² */
static inline LC3_FLOAT sqrf(LC3_FLOAT x) { return x * x; }

/* convenience wrappers around memmove */
static inline void move_float(LC3_FLOAT *dst, const LC3_FLOAT *src, LC3_INT len) {
  memmove(dst, src, len * sizeof(LC3_FLOAT));
}
static inline void move_int(LC3_INT *dst, const LC3_INT *src, LC3_INT len) {
  memmove(dst, src, len * sizeof(LC3_INT));
}

/* convenience wrappers around memset */
static inline void zero_float(LC3_FLOAT *x, LC3_INT len) {
  memset(x, 0, len * sizeof(LC3_FLOAT));
}
static inline void zero_int(LC3_INT *x, LC3_INT len) {
  memset(x, 0, len * sizeof(LC3_INT));
}

/* multiply float vectors element by element, in-place */
static inline void mult_vec(LC3_FLOAT *a, const LC3_FLOAT *b,
                            LC3_INT len) {
  LC3_INT i = 0;
  for (i = 0; i < len; i++) {
    a[i] *= b[i];
  }
}

/* multiply float vector with constant, in-place */
static inline void mult_const(LC3_FLOAT *a, LC3_FLOAT b, LC3_INT len) {
  LC3_INT i = 0;
  for (i = 0; i < len; i++) {
    a[i] *= b;
  }
}

/* sum of vector */
static inline LC3_FLOAT sum_vec(const LC3_FLOAT *x, LC3_INT len) {
  LC3_FLOAT sum = 0;
  LC3_INT i = 0;
  for (i = 0; i < len; i++) {
    sum += x[i];
  }
  return sum;
}

/* complex constructor */
static inline Complex cmplx(LC3_FLOAT r, LC3_FLOAT i) { return (Complex){r, i}; }

/* complex a + b */
static inline Complex cadd(Complex a, Complex b) {
  return cmplx(a.r + b.r, a.i + b.i);
}

/* complex a * b */
static inline Complex cmul(Complex a, Complex b) {
  return cmplx(a.r * b.r - a.i * b.i, a.i * b.r + a.r * b.i);
}

/* complex eᶦˣ */
static inline Complex cexpi(LC3_FLOAT x) { return cmplx(LC3_COS(x), LC3_SIN(x)); }

/* complex -x */
static inline Complex cneg(Complex x) { return cmplx(-x.r, -x.i); }

/* convert string to number. return true on success */
static inline bool str_to_int(const char *str, LC3_INT *value) {
  char *end = NULL;
  long v = str ? strtol(str, &end, 0) : 0;
  *value = (LC3_INT)v;
  return str && *end == 0 && v >= INT_MIN && v <= INT_MAX;
}

/* returns true if str ends with str ends with suffix. ignoring case. str may be
 * NULL */
static inline bool str_ends_with(const char *str, const char *suffix) {
  char *tmp = str ? strrchr(str, suffix[0]) : NULL;
  return tmp && !strcasecmp(tmp, suffix);
}

#endif
