/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Slatec's ZACON subroutine, originally translated from Fortran to
  C with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZACON_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZACON_IMPL_H

#include "f2c_d_sign_impl_.h"
#include "fmax_impl_.h"
#include "fortran_d1mach_impl_.h"
#include "min_impl_.h"
#include "slatec_zabs_impl_.h"
#include "slatec_zbinu_impl_.h"
#include "slatec_zbknu_impl_.h"
#include "slatec_zmlt_impl_.h"
#include "slatec_zs1s2_impl_.h"
#include <math.h> /* For math operations and constants */

/* zacon.f -- translated by f2c (version 20100827) */

/* DECK ZACON */
/* Subroutine */
static inline int slatec_zacon_impl_(double *zr, double *zi, double *fnu,
                                     int *kode, int *mr, int *n, double *yr,
                                     double *yi, int *nz, double *rl,
                                     double *fnul, double *tol, double *elim,
                                     double *alim) {
  /* Initialized data */
  int c_u_1 = 1;
  int c_u_2 = 2;
  double pi = 3.14159265358979324;
  double zeror = 0.;
  double coner = 1.;

  /* System generated locals */
  int i_u_1;

  /* Local variables */
  int i_u_;
  double fn;
  int nn, nw;
  double yy, c1i, c2i, c1m, as2, c1r, c2r, s1i, s2i, s1r, s2r, cki, arg, ckr,
      cpn;
  int iuf;
  double cyi[2], fmr, csr, azn, sgn;
  int inu;
  double bry[3], cyr[2], pti, spn, sti, zni, rzi, ptr, str, znr, rzr, sc1i,
      sc2i, sc1r, sc2r, cscl, cscr;
  double csrr[3], cssr[3], razn;
  int kflag;
  double ascle, bscle, csgni, csgnr, cspni, cspnr;

  /* ***BEGIN PROLOGUE  ZACON */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESH and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CACON-A, ZACON-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZACON APPLIES THE ANALYTIC CONTINUATION FORMULA */

  /*         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN) */
  /*                 MP=PI*MR*CMPLX(0.0,1.0) */

  /*     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT */
  /*     HALF Z PLANE */

  /* ***SEE ALSO  ZBESH, ZBESK */
  /* ***ROUTINES CALLED  D1MACH, ZABS, ZBINU, ZBKNU, ZMLT, ZS1S2 */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZACON */
  /*     COMPLEX CK,CONE,CSCL,CSCR,CSGN,CSPN,CY,CZERO,C1,C2,RZ,SC1,SC2,ST, */
  /*    *S1,S2,Y,Z,ZN */
  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZACON */
  *nz = 0;
  znr = -(*zr);
  zni = -(*zi);
  nn = *n;
  slatec_zbinu_impl_(&znr, &zni, fnu, kode, &nn, &yr[1], &yi[1], &nw, rl, fnul,
                     tol, elim, alim);
  if (nw < 0) {
    goto L90;
  }
  /* ----------------------------------------------------------------------- */
  /*     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION */
  /* ----------------------------------------------------------------------- */
  nn = min_impl_(2, *n);
  slatec_zbknu_impl_(&znr, &zni, fnu, kode, &nn, cyr, cyi, &nw, tol, elim,
                     alim);
  if (nw != 0) {
    goto L90;
  }
  s1r = cyr[0];
  s1i = cyi[0];
  fmr = (double)(*mr);
  sgn = -f2c_d_sign_impl_(&pi, &fmr);
  csgnr = zeror;
  csgni = sgn;
  if (*kode == 1) {
    goto L10;
  }
  yy = -zni;
  cpn = cos(yy);
  spn = sin(yy);
  slatec_zmlt_impl_(&csgnr, &csgni, &cpn, &spn, &csgnr, &csgni);
L10:
  /* ----------------------------------------------------------------------- */
  /*     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE */
  /*     WHEN FNU IS LARGE */
  /* ----------------------------------------------------------------------- */
  inu = (int)(*fnu);
  arg = (*fnu - inu) * sgn;
  cpn = cos(arg);
  spn = sin(arg);
  cspnr = cpn;
  cspni = spn;
  if (inu % 2 == 0) {
    goto L20;
  }
  cspnr = -cspnr;
  cspni = -cspni;
L20:
  iuf = 0;
  c1r = s1r;
  c1i = s1i;
  c2r = yr[1];
  c2i = yi[1];
  ascle = fortran_d1mach_impl_(&c_u_1) * 1e3 / *tol;
  if (*kode == 1) {
    goto L30;
  }
  slatec_zs1s2_impl_(&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, alim,
                     &iuf);
  *nz += nw;
  sc1r = c1r;
  sc1i = c1i;
L30:
  slatec_zmlt_impl_(&cspnr, &cspni, &c1r, &c1i, &str, &sti);
  slatec_zmlt_impl_(&csgnr, &csgni, &c2r, &c2i, &ptr, &pti);
  yr[1] = str + ptr;
  yi[1] = sti + pti;
  if (*n == 1) {
    return 0;
  }
  cspnr = -cspnr;
  cspni = -cspni;
  s2r = cyr[1];
  s2i = cyi[1];
  c1r = s2r;
  c1i = s2i;
  c2r = yr[2];
  c2i = yi[2];
  if (*kode == 1) {
    goto L40;
  }
  slatec_zs1s2_impl_(&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, alim,
                     &iuf);
  *nz += nw;
  sc2r = c1r;
  sc2i = c1i;
L40:
  slatec_zmlt_impl_(&cspnr, &cspni, &c1r, &c1i, &str, &sti);
  slatec_zmlt_impl_(&csgnr, &csgni, &c2r, &c2i, &ptr, &pti);
  yr[2] = str + ptr;
  yi[2] = sti + pti;
  if (*n == 2) {
    return 0;
  }
  cspnr = -cspnr;
  cspni = -cspni;
  azn = slatec_zabs_impl_(&znr, &zni);
  razn = 1. / azn;
  str = znr * razn;
  sti = -zni * razn;
  rzr = (str + str) * razn;
  rzi = (sti + sti) * razn;
  fn = *fnu + 1.;
  ckr = fn * rzr;
  cki = fn * rzi;
  /* ----------------------------------------------------------------------- */
  /*     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS */
  /* ----------------------------------------------------------------------- */
  cscl = 1. / *tol;
  cscr = *tol;
  cssr[0] = cscl;
  cssr[1] = coner;
  cssr[2] = cscr;
  csrr[0] = cscr;
  csrr[1] = coner;
  csrr[2] = cscl;
  bry[0] = ascle;
  bry[1] = 1. / ascle;
  bry[2] = fortran_d1mach_impl_(&c_u_2);
  as2 = slatec_zabs_impl_(&s2r, &s2i);
  kflag = 2;
  if (as2 > bry[0]) {
    goto L50;
  }
  kflag = 1;
  goto L60;
L50:
  if (as2 < bry[1]) {
    goto L60;
  }
  kflag = 3;
L60:
  bscle = bry[kflag - 1];
  s1r *= cssr[kflag - 1];
  s1i *= cssr[kflag - 1];
  s2r *= cssr[kflag - 1];
  s2i *= cssr[kflag - 1];
  csr = csrr[kflag - 1];
  i_u_1 = *n;
  for (i_u_ = 3; i_u_ <= i_u_1; ++i_u_) {
    str = s2r;
    sti = s2i;
    s2r = ckr * str - cki * sti + s1r;
    s2i = ckr * sti + cki * str + s1i;
    s1r = str;
    s1i = sti;
    c1r = s2r * csr;
    c1i = s2i * csr;
    str = c1r;
    sti = c1i;
    c2r = yr[i_u_];
    c2i = yi[i_u_];
    if (*kode == 1) {
      goto L70;
    }
    if (iuf < 0) {
      goto L70;
    }
    slatec_zs1s2_impl_(&znr, &zni, &c1r, &c1i, &c2r, &c2i, &nw, &ascle, alim,
                       &iuf);
    *nz += nw;
    sc1r = sc2r;
    sc1i = sc2i;
    sc2r = c1r;
    sc2i = c1i;
    if (iuf != 3) {
      goto L70;
    }
    iuf = -4;
    s1r = sc1r * cssr[kflag - 1];
    s1i = sc1i * cssr[kflag - 1];
    s2r = sc2r * cssr[kflag - 1];
    s2i = sc2i * cssr[kflag - 1];
    str = sc2r;
    sti = sc2i;
  L70:
    ptr = cspnr * c1r - cspni * c1i;
    pti = cspnr * c1i + cspni * c1r;
    yr[i_u_] = ptr + csgnr * c2r - csgni * c2i;
    yi[i_u_] = pti + csgnr * c2i + csgni * c2r;
    ckr += rzr;
    cki += rzi;
    cspnr = -cspnr;
    cspni = -cspni;
    if (kflag >= 3) {
      goto L80;
    }
    ptr = fabs(c1r);
    pti = fabs(c1i);
    c1m = fmax_impl_(ptr, pti);
    if (c1m <= bscle) {
      goto L80;
    }
    ++kflag;
    bscle = bry[kflag - 1];
    s1r *= csr;
    s1i *= csr;
    s2r = str;
    s2i = sti;
    s1r *= cssr[kflag - 1];
    s1i *= cssr[kflag - 1];
    s2r *= cssr[kflag - 1];
    s2i *= cssr[kflag - 1];
    csr = csrr[kflag - 1];
  L80:;
  }
  return 0;
L90:
  *nz = -1;
  if (nw == -2) {
    *nz = -2;
  }
  return 0;
} /* slatec_zacon_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZACON_IMPL_H */