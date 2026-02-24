/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Slatec's ZUCHK subroutine, originally translated from Fortran to
  C with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZUCHK_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZUCHK_IMPL_H

#include "fmax_impl_.h"
#include "fmin_impl_.h"
#include <math.h> /* For math operations and constants */

/* zuchk.f -- translated by f2c (version 20100827) */

/* DECK ZUCHK */
/* Subroutine */
static inline int slatec_zuchk_impl_(double *yr, double *yi, int *nz,
                                     double *ascle, double *tol) {
  double wi, ss, st, wr;

  /* ***BEGIN PROLOGUE  ZUCHK */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to SERI, ZUOIK, ZUNK1, ZUNK2, ZUNI1, ZUNI2 and */
  /*            ZKSCL */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CUCHK-A, ZUCHK-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN */
  /*      EXP(-ALIM)=ASCLE=1.0E+3*D1MACH(1)/TOL. THE TEST IS MADE TO SEE */
  /*      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW */
  /*      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED */
  /*      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE */
  /*      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE */
  /*      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED. */

  /* ***SEE ALSO  SERI, ZKSCL, ZUNI1, ZUNI2, ZUNK1, ZUNK2, ZUOIK */
  /* ***ROUTINES CALLED  (NONE) */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   ??????  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZUCHK */

  /*     COMPLEX Y */
  /* ***FIRST EXECUTABLE STATEMENT  ZUCHK */
  *nz = 0;
  wr = fabs(*yr);
  wi = fabs(*yi);
  st = fmin_impl_(wr, wi);
  if (st > *ascle) {
    return 0;
  }
  ss = fmax_impl_(wr, wi);
  st /= *tol;
  if (ss < st) {
    *nz = 1;
  }
  return 0;
} /* slatec_zuchk_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZUCHK_IMPL_H */