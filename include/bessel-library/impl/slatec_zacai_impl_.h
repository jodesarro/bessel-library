/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Slatec's ZACAI subroutine, originally translated from Fortran to
  C with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZACAI_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZACAI_IMPL_H

#include "f2c_d_sign_impl_.h"
#include "fortran_d1mach_impl_.h"
#include "slatec_zabs_impl_.h"
#include "slatec_zasyi_impl_.h"
#include "slatec_zbknu_impl_.h"
#include "slatec_zmlri_impl_.h"
#include "slatec_zs1s2_impl_.h"
#include "slatec_zseri_impl_.h"
#include <math.h> /* For math operations and constants */

/* zacai.f -- translated by f2c (version 20100827) */

/* DECK ZACAI */
/* Subroutine */
static inline int slatec_zacai_impl_(double *zr, double *zi, double *fnu,
                                     int *kode, int *mr, int *n, double *yr,
                                     double *yi, int *nz, double *rl,
                                     double *tol, double *elim, double *alim) {
  /* Initialized data */
  static int c_u_1 = 1;
  static double pi = 3.14159265358979324;

  /* Local variables */
  static double az;
  static int nn, nw;
  static double yy, c1i, c2i, c1r, c2r, arg;
  static int iuf;
  static double cyi[2], fmr, sgn;
  static int inu;
  static double cyr[2], zni, znr, dfnu;
  static double ascle, csgni, csgnr, cspni, cspnr;

  /* ***BEGIN PROLOGUE  ZACAI */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZAIRY */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CACAI-A, ZACAI-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA */

  /*         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN) */
  /*                 MP=PI*MR*CMPLX(0.0,1.0) */

  /*     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT */
  /*     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1. */
  /*     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND */
  /*     RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT IF ZACON */
  /*     IS CALLED FROM ZAIRY. */

  /* ***SEE ALSO  ZAIRY */
  /* ***ROUTINES CALLED  D1MACH, ZABS, ZASYI, ZBKNU, ZMLRI, ZS1S2, ZSERI */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZACAI */
  /*     COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY */
  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZACAI */
  *nz = 0;
  znr = -(*zr);
  zni = -(*zi);
  az = slatec_zabs_impl_(zr, zi);
  nn = *n;
  dfnu = *fnu + (*n - 1);
  if (az <= 2.) {
    goto L10;
  }
  if (az * az * .25 > dfnu + 1.) {
    goto L20;
  }
L10:
  /* ----------------------------------------------------------------------- */
  /*     POWER SERIES FOR THE I FUNCTION */
  /* ----------------------------------------------------------------------- */
  slatec_zseri_impl_(&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, tol, elim,
                     alim);
  goto L40;
L20:
  if (az < *rl) {
    goto L30;
  }
  /* ----------------------------------------------------------------------- */
  /*     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION */
  /* ----------------------------------------------------------------------- */
  slatec_zasyi_impl_(&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, rl, tol,
                     elim, alim);
  if (nw < 0) {
    goto L80;
  }
  goto L40;
L30:
  /* ----------------------------------------------------------------------- */
  /*     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION */
  /* ----------------------------------------------------------------------- */
  slatec_zmlri_impl_(&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, tol);
  if (nw < 0) {
    goto L80;
  }
L40:
  /* ----------------------------------------------------------------------- */
  /*     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION */
  /* ----------------------------------------------------------------------- */
  slatec_zbknu_impl_(&znr, &zni, fnu, kode, &c_u_1, cyr, cyi, &nw, tol, elim,
                     alim);
  if (nw != 0) {
    goto L80;
  }
  fmr = (double)(*mr);
  sgn = -f2c_d_sign_impl_(&pi, &fmr);
  csgnr = 0.;
  csgni = sgn;
  if (*kode == 1) {
    goto L50;
  }
  yy = -zni;
  csgnr = -csgni * sin(yy);
  csgni *= cos(yy);
L50:
  /* ----------------------------------------------------------------------- */
  /*     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE */
  /*     WHEN FNU IS LARGE */
  /* ----------------------------------------------------------------------- */
  inu = (int)(*fnu);
  arg = (*fnu - inu) * sgn;
  cspnr = cos(arg);
  cspni = sin(arg);
  if (inu % 2 == 0) {
    goto L60;
  }
  cspnr = -cspnr;
  cspni = -cspni;
L60:
  c1r = cyr[0];
  c1i = cyi[0];
  c2r = yr[1];
  c2i = yi[1];
  if (*kode == 1) {
    goto L70;
  }
  iuf = 0;
  ascle = fortran_d1mach_impl_(&c_u_1) * 1e3 / *tol;
  slatec_zs1s2_impl_(&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, alim,
                     &iuf);
  *nz += nw;
L70:
  yr[1] = cspnr * c1r - cspni * c1i + csgnr * c2r - csgni * c2i;
  yi[1] = cspnr * c1i + cspni * c1r + csgnr * c2i + csgni * c2r;
  return 0;
L80:
  *nz = -1;
  if (nw == -2) {
    *nz = -2;
  }
  return 0;
} /* slatec_zacai_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZACAI_IMPL_H */