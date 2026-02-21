/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Slatec's ZLOG subroutine, originally translated from Fortran to C
  with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZLOG_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZLOG_IMPL_H

#include "slatec_zabs_impl_.h"
#include <math.h> /* For math operations and constants */

/* zlog.f -- translated by f2c (version 20100827) */

/* DECK ZLOG */
/* Subroutine */
static inline int slatec_zlog_impl_(double *ar, double *ai, double *br,
                                    double *bi, int *ierr) {
  /* Initialized data */
  static double dpi = 3.141592653589793238462643383;
  static double dhpi = 1.570796326794896619231321696;

  /* Local variables */
  static double zm;
  static double dtheta;

  /* ***BEGIN PROLOGUE  ZLOG */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and */
  /*            ZBIRY */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (ZLOG-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     DOUBLE PRECISION COMPLEX LOGARITHM B=CLOG(A) */
  /*     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0) */
  /* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY */
  /* ***ROUTINES CALLED  ZABS */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZLOG */
  /* ***FIRST EXECUTABLE STATEMENT  ZLOG */
  *ierr = 0;
  if (*ar == 0.) {
    goto L10;
  }
  if (*ai == 0.) {
    goto L20;
  }
  dtheta = atan(*ai / *ar);
  if (dtheta <= 0.) {
    goto L40;
  }
  if (*ar < 0.) {
    dtheta -= dpi;
  }
  goto L50;
L10:
  if (*ai == 0.) {
    goto L60;
  }
  *bi = dhpi;
  *br = log((fabs(*ai)));
  if (*ai < 0.) {
    *bi = -(*bi);
  }
  return 0;
L20:
  if (*ar > 0.) {
    goto L30;
  }
  *br = log((fabs(*ar)));
  *bi = dpi;
  return 0;
L30:
  *br = log(*ar);
  *bi = 0.;
  return 0;
L40:
  if (*ar < 0.) {
    dtheta += dpi;
  }
L50:
  zm = slatec_zabs_impl_(ar, ai);
  *br = log(zm);
  *bi = dtheta;
  return 0;
L60:
  *ierr = 1;
  return 0;
} /* slatec_zlog_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZLOG_IMPL_H */