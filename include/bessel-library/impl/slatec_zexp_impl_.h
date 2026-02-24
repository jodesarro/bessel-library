/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Slatec's ZEXP subroutine, originally translated from Fortran to C
  with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZEXP_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZEXP_IMPL_H

#include <math.h> /* For math operations and constants */

/* zexp.f -- translated by f2c (version 20100827) */

/* DECK ZEXP */
/* Subroutine */
static inline int slatec_zexp_impl_(double *ar, double *ai, double *br,
                                    double *bi) {
  /* Local variables */
  double ca, cb, zm;

  /* ***BEGIN PROLOGUE  ZEXP */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and */
  /*            ZBIRY */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (ZEXP-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B=EXP(A) */

  /* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY */
  /* ***ROUTINES CALLED  (NONE) */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZEXP */
  /* ***FIRST EXECUTABLE STATEMENT  ZEXP */
  zm = exp(*ar);
  ca = zm * cos(*ai);
  cb = zm * sin(*ai);
  *br = ca;
  *bi = cb;
  return 0;
} /* slatec_zexp_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZEXP_IMPL_H */