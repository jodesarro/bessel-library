/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Slatec's ZMLT subroutine, originally translated from Fortran to C
  with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZMLT_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZMLT_IMPL_H

/* zmlt.f -- translated by f2c (version 20100827) */

/* DECK ZMLT */
/* Subroutine */
static inline int slatec_zmlt_impl_(double *ar, double *ai, double *br,
                                    double *bi, double *cr, double *ci) {
  double ca, cb;

  /* ***BEGIN PROLOGUE  ZMLT */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and */
  /*            ZBIRY */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (ZMLT-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B. */

  /* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY */
  /* ***ROUTINES CALLED  (NONE) */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZMLT */
  /* ***FIRST EXECUTABLE STATEMENT  ZMLT */
  ca = *ar * *br - *ai * *bi;
  cb = *ar * *bi + *ai * *br;
  *cr = ca;
  *ci = cb;
  return 0;
} /* slatec_zmlt_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZMLT_IMPL_H */