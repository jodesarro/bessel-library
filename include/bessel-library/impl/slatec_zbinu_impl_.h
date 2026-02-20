/*
  Bessel Library: A C library with routines for computing Bessel functions

  File: include/bessel-library/impl/slatec_zbinu_impl_.h
  Language standards: C99
  References: include/bessel-library/references.txt
  License: include/bessel-library/license.txt
  Repository: <https://github.com/jodesarro/bessel-library>

  Description: Slatec's ZBINU subroutine, originally translated from Fortran to
  C with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZBINU_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZBINU_IMPL_H

#include "max_impl_.h"
#include "slatec_zabs_impl_.h"
#include "slatec_zasyi_impl_.h"
#include "slatec_zbuni_impl_.h"
#include "slatec_zmlri_impl_.h"
#include "slatec_zseri_impl_.h"
#include "slatec_zuoik_impl_.h"
#include "slatec_zwrsk_impl_.h"
#include <stdlib.h> /* For abs() */

/* zbinu.f -- translated by f2c (version 20100827) */

/* DECK ZBINU */
/* Subroutine */
static inline int slatec_zbinu_impl_(double *zr, double *zi, double *fnu,
                                     int *kode, int *n, double *cyr,
                                     double *cyi, int *nz, double *rl,
                                     double *fnul, double *tol, double *elim,
                                     double *alim) {
  /* Initialized data */
  static int c_u_1 = 1;
  static int c_u_2 = 2;
  static double zeror = 0.;
  static double zeroi = 0.;

  /* System generated locals */
  int i_u_1;

  /* Local variables */
  static int i_u_;
  static double az;
  static int nn, nw;
  static double cwi[2], cwr[2];
  static int nui, inw;
  static double dfnu;
  static int nlast;

  /* ***BEGIN PROLOGUE  ZBINU */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK and ZBIRY */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CBINU-A, ZBINU-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE */

  /* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBIRY */
  /* ***ROUTINES CALLED  ZABS, ZASYI, ZBUNI, ZMLRI, ZSERI, ZUOIK, ZWRSK */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZBINU */
  /* Parameter adjustments */
  --cyi;
  --cyr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZBINU */
  *nz = 0;
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
  /*     POWER SERIES */
  /* ----------------------------------------------------------------------- */
  slatec_zseri_impl_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, tol, elim,
                     alim);
  inw = abs(nw);
  *nz += inw;
  nn -= inw;
  if (nn == 0) {
    return 0;
  }
  if (nw >= 0) {
    goto L120;
  }
  dfnu = *fnu + (nn - 1);
L20:
  if (az < *rl) {
    goto L40;
  }
  if (dfnu <= 1.) {
    goto L30;
  }
  if (az + az < dfnu * dfnu) {
    goto L50;
  }
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR LARGE Z */
/* ----------------------------------------------------------------------- */
L30:
  slatec_zasyi_impl_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, rl, tol,
                     elim, alim);
  if (nw < 0) {
    goto L130;
  }
  goto L120;
L40:
  if (dfnu <= 1.) {
    goto L70;
  }
L50:
  /* ----------------------------------------------------------------------- */
  /*     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM */
  /* ----------------------------------------------------------------------- */
  slatec_zuoik_impl_(zr, zi, fnu, kode, &c_u_1, &nn, &cyr[1], &cyi[1], &nw, tol,
                     elim, alim);
  if (nw < 0) {
    goto L130;
  }
  *nz += nw;
  nn -= nw;
  if (nn == 0) {
    return 0;
  }
  dfnu = *fnu + (nn - 1);
  if (dfnu > *fnul) {
    goto L110;
  }
  if (az > *fnul) {
    goto L110;
  }
L60:
  if (az > *rl) {
    goto L80;
  }
L70:
  /* ----------------------------------------------------------------------- */
  /*     MILLER ALGORITHM NORMALIZED BY THE SERIES */
  /* ----------------------------------------------------------------------- */
  slatec_zmlri_impl_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, tol);
  if (nw < 0) {
    goto L130;
  }
  goto L120;
L80:
  /* ----------------------------------------------------------------------- */
  /*     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN */
  /* ----------------------------------------------------------------------- */
  /* ----------------------------------------------------------------------- */
  /*     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN */
  /* ----------------------------------------------------------------------- */
  slatec_zuoik_impl_(zr, zi, fnu, kode, &c_u_2, &c_u_2, cwr, cwi, &nw, tol,
                     elim, alim);
  if (nw >= 0) {
    goto L100;
  }
  *nz = nn;
  i_u_1 = nn;
  for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
    cyr[i_u_] = zeror;
    cyi[i_u_] = zeroi;
    /* L90: */
  }
  return 0;
L100:
  if (nw > 0) {
    goto L130;
  }
  slatec_zwrsk_impl_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, cwr, cwi,
                     tol, elim, alim);
  if (nw < 0) {
    goto L130;
  }
  goto L120;
L110:
  /* ----------------------------------------------------------------------- */
  /*     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD */
  /* ----------------------------------------------------------------------- */
  nui = (int)(*fnul - dfnu + 1);
  nui = max_impl_(nui, 0);
  slatec_zbuni_impl_(zr, zi, fnu, kode, &nn, &cyr[1], &cyi[1], &nw, &nui,
                     &nlast, fnul, tol, elim, alim);
  if (nw < 0) {
    goto L130;
  }
  *nz += nw;
  if (nlast == 0) {
    goto L120;
  }
  nn = nlast;
  goto L60;
L120:
  return 0;
L130:
  *nz = -1;
  if (nw == -2) {
    *nz = -2;
  }
  return 0;
} /* slatec_zbinu_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZBINU_IMPL_H */