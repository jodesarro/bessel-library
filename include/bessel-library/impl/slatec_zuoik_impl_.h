/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Slatec's ZUOIK subroutine, originally translated from Fortran to
  C with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZUOIK_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZUOIK_IMPL_H

#include "fmax_impl_.h"
#include "fortran_d1mach_impl_.h"
#include "slatec_zabs_impl_.h"
#include "slatec_zlog_impl_.h"
#include "slatec_zuchk_impl_.h"
#include "slatec_zunhj_impl_.h"
#include "slatec_zunik_impl_.h"
#include <math.h> /* For math operations and constants */

/* zuoik.f -- translated by f2c (version 20100827) */

/* DECK ZUOIK */
/* Subroutine */
static inline int slatec_zuoik_impl_(double *zr, double *zi, double *fnu,
                                     int *kode, int *ikflg, int *n, double *yr,
                                     double *yi, int *nuf, double *tol,
                                     double *elim, double *alim) {
  /* Initialized data */
  int c_u_1 = 1;
  double zeror = 0.;
  double zeroi = 0.;
  double aic = 1.265512123484645396;

  /* System generated locals */
  int i_u_1;

  /* Local variables */
  int i_u_;
  double ax, ay;
  int nn, nw;
  double fnn, gnn, zbi, czi, gnu, zbr, czr, rcz, sti, zni, zri, str, znr, zrr,
      aarg, aphi, argi, phii, argr;
  int idum;
  double phir;
  int init;
  double sumi, sumr, ascle;
  int iform;
  double asumi, bsumi, cwrki[16];
  double asumr, bsumr, cwrkr[16];
  double zeta1i, zeta2i, zeta1r, zeta2r;

  /* ***BEGIN PROLOGUE  ZUOIK */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESH, ZBESI and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CUOIK-A, ZUOIK-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC */
  /*     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM */
  /*     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW */
  /*     WHERE ALIM.LT.ELIM. IF THE MAGNITUDE, BASED ON THE LEADING */
  /*     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN */
  /*     THE RESULT IS ON SCALE. IF NOT, THEN A REFINED TEST USING OTHER */
  /*     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE */
  /*     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)= */
  /*     EXP(-ELIM)/TOL */

  /*     IKFLG=1 MEANS THE I SEQUENCE IS TESTED */
  /*          =2 MEANS THE K SEQUENCE IS TESTED */
  /*     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE */
  /*         =-1 MEANS AN OVERFLOW WOULD OCCUR */
  /*     IKFLG=1 AND NUF.GT.0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO */
  /*             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE */
  /*     IKFLG=2 AND NUF.EQ.N MEANS ALL Y VALUES WERE SET TO ZERO */
  /*     IKFLG=2 AND 0.LT.NUF.LT.N NOT CONSIDERED. Y MUST BE SET BY */
  /*             ANOTHER ROUTINE */

  /* ***SEE ALSO  ZBESH, ZBESI, ZBESK */
  /* ***ROUTINES CALLED  D1MACH, ZABS, ZLOG, ZUCHK, ZUNHJ, ZUNIK */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /*   930122  Added ZLOG to EXTERNAL statement.  (RWC) */
  /* ***END PROLOGUE  ZUOIK */
  /*     COMPLEX ARG,ASUM,BSUM,CWRK,CZ,CZERO,PHI,SUM,Y,Z,ZB,ZETA1,ZETA2,ZN, */
  /*    *ZR */
  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZUOIK */
  *nuf = 0;
  nn = *n;
  zrr = *zr;
  zri = *zi;
  if (*zr >= 0.) {
    goto L10;
  }
  zrr = -(*zr);
  zri = -(*zi);
L10:
  zbr = zrr;
  zbi = zri;
  ax = fabs(*zr) * 1.7321;
  ay = fabs(*zi);
  iform = 1;
  if (ay > ax) {
    iform = 2;
  }
  gnu = fmax_impl_(*fnu, 1.);
  if (*ikflg == 1) {
    goto L20;
  }
  fnn = (double)nn;
  gnn = *fnu + fnn - 1.;
  gnu = fmax_impl_(gnn, fnn);
L20:
  /* ----------------------------------------------------------------------- */
  /*     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE */
  /*     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET */
  /*     THE SIGN OF THE IMAGINARY PART CORRECT. */
  /* ----------------------------------------------------------------------- */
  if (iform == 2) {
    goto L30;
  }
  init = 0;
  slatec_zunik_impl_(&zrr, &zri, &gnu, ikflg, &c_u_1, tol, &init, &phir, &phii,
                     &zeta1r, &zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr,
                     cwrki);
  czr = -zeta1r + zeta2r;
  czi = -zeta1i + zeta2i;
  goto L50;
L30:
  znr = zri;
  zni = -zrr;
  if (*zi > 0.) {
    goto L40;
  }
  znr = -znr;
L40:
  slatec_zunhj_impl_(&znr, &zni, &gnu, &c_u_1, tol, &phir, &phii, &argr, &argi,
                     &zeta1r, &zeta1i, &zeta2r, &zeta2i, &asumr, &asumi, &bsumr,
                     &bsumi);
  czr = -zeta1r + zeta2r;
  czi = -zeta1i + zeta2i;
  aarg = slatec_zabs_impl_(&argr, &argi);
L50:
  if (*kode == 1) {
    goto L60;
  }
  czr -= zbr;
  czi -= zbi;
L60:
  if (*ikflg == 1) {
    goto L70;
  }
  czr = -czr;
  czi = -czi;
L70:
  aphi = slatec_zabs_impl_(&phir, &phii);
  rcz = czr;
  /* ----------------------------------------------------------------------- */
  /*     OVERFLOW TEST */
  /* ----------------------------------------------------------------------- */
  if (rcz > *elim) {
    goto L210;
  }
  if (rcz < *alim) {
    goto L80;
  }
  rcz += log(aphi);
  if (iform == 2) {
    rcz = rcz - log(aarg) * .25 - aic;
  }
  if (rcz > *elim) {
    goto L210;
  }
  goto L130;
L80:
  /* ----------------------------------------------------------------------- */
  /*     UNDERFLOW TEST */
  /* ----------------------------------------------------------------------- */
  if (rcz < -(*elim)) {
    goto L90;
  }
  if (rcz > -(*alim)) {
    goto L130;
  }
  rcz += log(aphi);
  if (iform == 2) {
    rcz = rcz - log(aarg) * .25 - aic;
  }
  if (rcz > -(*elim)) {
    goto L110;
  }
L90:
  i_u_1 = nn;
  for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
    yr[i_u_] = zeror;
    yi[i_u_] = zeroi;
    /* L100: */
  }
  *nuf = nn;
  return 0;
L110:
  ascle = fortran_d1mach_impl_(&c_u_1) * 1e3 / *tol;
  slatec_zlog_impl_(&phir, &phii, &str, &sti, &idum);
  czr += str;
  czi += sti;
  if (iform == 1) {
    goto L120;
  }
  slatec_zlog_impl_(&argr, &argi, &str, &sti, &idum);
  czr = czr - str * .25 - aic;
  czi -= sti * .25;
L120:
  ax = exp(rcz) / *tol;
  ay = czi;
  czr = ax * cos(ay);
  czi = ax * sin(ay);
  slatec_zuchk_impl_(&czr, &czi, &nw, &ascle, tol);
  if (nw != 0) {
    goto L90;
  }
L130:
  if (*ikflg == 2) {
    return 0;
  }
  if (*n == 1) {
    return 0;
  }
/* ----------------------------------------------------------------------- */
/*     SET UNDERFLOWS ON I SEQUENCE */
/* ----------------------------------------------------------------------- */
L140:
  gnu = *fnu + (nn - 1);
  if (iform == 2) {
    goto L150;
  }
  init = 0;
  slatec_zunik_impl_(&zrr, &zri, &gnu, ikflg, &c_u_1, tol, &init, &phir, &phii,
                     &zeta1r, &zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr,
                     cwrki);
  czr = -zeta1r + zeta2r;
  czi = -zeta1i + zeta2i;
  goto L160;
L150:
  slatec_zunhj_impl_(&znr, &zni, &gnu, &c_u_1, tol, &phir, &phii, &argr, &argi,
                     &zeta1r, &zeta1i, &zeta2r, &zeta2i, &asumr, &asumi, &bsumr,
                     &bsumi);
  czr = -zeta1r + zeta2r;
  czi = -zeta1i + zeta2i;
  aarg = slatec_zabs_impl_(&argr, &argi);
L160:
  if (*kode == 1) {
    goto L170;
  }
  czr -= zbr;
  czi -= zbi;
L170:
  aphi = slatec_zabs_impl_(&phir, &phii);
  rcz = czr;
  if (rcz < -(*elim)) {
    goto L180;
  }
  if (rcz > -(*alim)) {
    return 0;
  }
  rcz += log(aphi);
  if (iform == 2) {
    rcz = rcz - log(aarg) * .25 - aic;
  }
  if (rcz > -(*elim)) {
    goto L190;
  }
L180:
  yr[nn] = zeror;
  yi[nn] = zeroi;
  --nn;
  ++(*nuf);
  if (nn == 0) {
    return 0;
  }
  goto L140;
L190:
  ascle = fortran_d1mach_impl_(&c_u_1) * 1e3 / *tol;
  slatec_zlog_impl_(&phir, &phii, &str, &sti, &idum);
  czr += str;
  czi += sti;
  if (iform == 1) {
    goto L200;
  }
  slatec_zlog_impl_(&argr, &argi, &str, &sti, &idum);
  czr = czr - str * .25 - aic;
  czi -= sti * .25;
L200:
  ax = exp(rcz) / *tol;
  ay = czi;
  czr = ax * cos(ay);
  czi = ax * sin(ay);
  slatec_zuchk_impl_(&czr, &czi, &nw, &ascle, tol);
  if (nw != 0) {
    goto L180;
  }
  return 0;
L210:
  *nuf = -1;
  return 0;
} /* slatec_zuoik_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZUOIK_IMPL_H */