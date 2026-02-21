/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Slatec's ZUNK1 subroutine, originally translated from Fortran to
  C with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZUNK1_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZUNK1_IMPL_H

#include "f2c_d_sign_impl_.h"
#include "fmax_impl_.h"
#include "fortran_d1mach_impl_.h"
#include "slatec_zabs_impl_.h"
#include "slatec_zs1s2_impl_.h"
#include "slatec_zuchk_impl_.h"
#include "slatec_zunik_impl_.h"
#include <math.h> /* For math operations and constants */

/* zunk1.f -- translated by f2c (version 20100827) */

/* DECK ZUNK1 */
/* Subroutine */
static inline int slatec_zunk1_impl_(double *zr, double *zi, double *fnu,
                                     int *kode, int *mr, int *n, double *yr,
                                     double *yi, int *nz, double *tol,
                                     double *elim, double *alim) {
  /* Initialized data */
  static int c_u_0 = 0;
  static int c_u_1 = 1;
  static int c_u_2 = 2;
  static double zeror = 0.;
  static double zeroi = 0.;
  static double coner = 1.;
  static double pi = 3.14159265358979324;

  /* System generated locals */
  int i_u_1;

  /* Local variables */
  static int i_u_, j, k, m, ib, ic;
  static double fn;
  static int il, kk, nw;
  static double c1i, c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r, ang, asc, cki,
      fnf;
  static int ifn;
  static double ckr;
  static int iuf;
  static double cyi[2], fmr, csr, sgn;
  static int inu;
  static double bry[3], cyr[2], sti, rzi, zri, str, rzr, zrr, aphi, cscl,
      phii[2], crsc;
  static double phir[2];
  static int init[2];
  static double csrr[3], cssr[3], rast, sumi[2], razr;
  static double sumr[2];
  static int iflag, kflag;
  static double ascle;
  static int kdflg;
  static double phidi;
  static int ipard;
  static double csgni, phidr;
  static int initd;
  static double cspni, cwrki[48] /* was [16][3] */, sumdi;
  static double cspnr, cwrkr[48] /* was [16][3] */, sumdr;
  static double zeta1i[2], zeta2i[2], zet1di, zet2di, zeta1r[2], zeta2r[2],
      zet1dr, zet2dr;

  /* ***BEGIN PROLOGUE  ZUNK1 */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CUNK1-A, ZUNK1-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE */
  /*     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE */
  /*     UNIFORM ASYMPTOTIC EXPANSION. */
  /*     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION. */
  /*     NZ=-1 MEANS AN OVERFLOW WILL OCCUR */

  /* ***SEE ALSO  ZBESK */
  /* ***ROUTINES CALLED  D1MACH, ZABS, ZS1S2, ZUCHK, ZUNIK */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZUNK1 */
  /*     COMPLEX CFN,CK,CONE,CRSC,CS,CSCL,CSGN,CSPN,CSR,CSS,CWRK,CY,CZERO, */
  /*    *C1,C2,PHI,PHID,RZ,SUM,SUMD,S1,S2,Y,Z,ZETA1,ZETA1D,ZETA2,ZETA2D,ZR */
  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZUNK1 */
  kdflg = 1;
  *nz = 0;
  /* ----------------------------------------------------------------------- */
  /*     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN */
  /*     THE UNDERFLOW LIMIT */
  /* ----------------------------------------------------------------------- */
  cscl = 1. / *tol;
  crsc = *tol;
  cssr[0] = cscl;
  cssr[1] = coner;
  cssr[2] = crsc;
  csrr[0] = crsc;
  csrr[1] = coner;
  csrr[2] = cscl;
  bry[0] = fortran_d1mach_impl_(&c_u_1) * 1e3 / *tol;
  bry[1] = 1. / bry[0];
  bry[2] = fortran_d1mach_impl_(&c_u_2);
  zrr = *zr;
  zri = *zi;
  if (*zr >= 0.) {
    goto L10;
  }
  zrr = -(*zr);
  zri = -(*zi);
L10:
  j = 2;
  i_u_1 = *n;
  for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
    /* -----------------------------------------------------------------------
     */
    /*     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J */
    /* -----------------------------------------------------------------------
     */
    j = 3 - j;
    fn = *fnu + (i_u_ - 1);
    init[j - 1] = 0;
    slatec_zunik_impl_(&zrr, &zri, &fn, &c_u_2, &c_u_0, tol, &init[j - 1],
                       &phir[j - 1], &phii[j - 1], &zeta1r[j - 1],
                       &zeta1i[j - 1], &zeta2r[j - 1], &zeta2i[j - 1],
                       &sumr[j - 1], &sumi[j - 1], &cwrkr[(j << 4) - 16],
                       &cwrki[(j << 4) - 16]);
    if (*kode == 1) {
      goto L20;
    }
    str = zrr + zeta2r[j - 1];
    sti = zri + zeta2i[j - 1];
    rast = fn / slatec_zabs_impl_(&str, &sti);
    str = str * rast * rast;
    sti = -sti * rast * rast;
    s1r = zeta1r[j - 1] - str;
    s1i = zeta1i[j - 1] - sti;
    goto L30;
  L20:
    s1r = zeta1r[j - 1] - zeta2r[j - 1];
    s1i = zeta1i[j - 1] - zeta2i[j - 1];
  L30:
    rs1 = s1r;
    /* -----------------------------------------------------------------------
     */
    /*     TEST FOR UNDERFLOW AND OVERFLOW */
    /* -----------------------------------------------------------------------
     */
    if (fabs(rs1) > *elim) {
      goto L60;
    }
    if (kdflg == 1) {
      kflag = 2;
    }
    if (fabs(rs1) < *alim) {
      goto L40;
    }
    /* -----------------------------------------------------------------------
     */
    /*     REFINE  TEST AND SCALE */
    /* -----------------------------------------------------------------------
     */
    aphi = slatec_zabs_impl_(&phir[j - 1], &phii[j - 1]);
    rs1 += log(aphi);
    if (fabs(rs1) > *elim) {
      goto L60;
    }
    if (kdflg == 1) {
      kflag = 1;
    }
    if (rs1 < 0.) {
      goto L40;
    }
    if (kdflg == 1) {
      kflag = 3;
    }
  L40:
    /* -----------------------------------------------------------------------
     */
    /*     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR */
    /*     EXPONENT EXTREMES */
    /* -----------------------------------------------------------------------
     */
    s2r = phir[j - 1] * sumr[j - 1] - phii[j - 1] * sumi[j - 1];
    s2i = phir[j - 1] * sumi[j - 1] + phii[j - 1] * sumr[j - 1];
    str = exp(s1r) * cssr[kflag - 1];
    s1r = str * cos(s1i);
    s1i = str * sin(s1i);
    str = s2r * s1r - s2i * s1i;
    s2i = s1r * s2i + s2r * s1i;
    s2r = str;
    if (kflag != 1) {
      goto L50;
    }
    slatec_zuchk_impl_(&s2r, &s2i, &nw, bry, tol);
    if (nw != 0) {
      goto L60;
    }
  L50:
    cyr[kdflg - 1] = s2r;
    cyi[kdflg - 1] = s2i;
    yr[i_u_] = s2r * csrr[kflag - 1];
    yi[i_u_] = s2i * csrr[kflag - 1];
    if (kdflg == 2) {
      goto L75;
    }
    kdflg = 2;
    goto L70;
  L60:
    if (rs1 > 0.) {
      goto L300;
    }
    /* -----------------------------------------------------------------------
     */
    /*     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW */
    /* -----------------------------------------------------------------------
     */
    if (*zr < 0.) {
      goto L300;
    }
    kdflg = 1;
    yr[i_u_] = zeror;
    yi[i_u_] = zeroi;
    ++(*nz);
    if (i_u_ == 1) {
      goto L70;
    }
    if (yr[i_u_ - 1] == zeror && yi[i_u_ - 1] == zeroi) {
      goto L70;
    }
    yr[i_u_ - 1] = zeror;
    yi[i_u_ - 1] = zeroi;
    ++(*nz);
  L70:;
  }
  i_u_ = *n;
L75:
  razr = 1. / slatec_zabs_impl_(&zrr, &zri);
  str = zrr * razr;
  sti = -zri * razr;
  rzr = (str + str) * razr;
  rzi = (sti + sti) * razr;
  ckr = fn * rzr;
  cki = fn * rzi;
  ib = i_u_ + 1;
  if (*n < ib) {
    goto L160;
  }
  /* ----------------------------------------------------------------------- */
  /*     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO */
  /*     ON UNDERFLOW. */
  /* ----------------------------------------------------------------------- */
  fn = *fnu + (*n - 1);
  ipard = 1;
  if (*mr != 0) {
    ipard = 0;
  }
  initd = 0;
  slatec_zunik_impl_(&zrr, &zri, &fn, &c_u_2, &ipard, tol, &initd, &phidr,
                     &phidi, &zet1dr, &zet1di, &zet2dr, &zet2di, &sumdr, &sumdi,
                     &cwrkr[32], &cwrki[32]);
  if (*kode == 1) {
    goto L80;
  }
  str = zrr + zet2dr;
  sti = zri + zet2di;
  rast = fn / slatec_zabs_impl_(&str, &sti);
  str = str * rast * rast;
  sti = -sti * rast * rast;
  s1r = zet1dr - str;
  s1i = zet1di - sti;
  goto L90;
L80:
  s1r = zet1dr - zet2dr;
  s1i = zet1di - zet2di;
L90:
  rs1 = s1r;
  if (fabs(rs1) > *elim) {
    goto L95;
  }
  if (fabs(rs1) < *alim) {
    goto L100;
  }
  /* ----------------------------------------------------------------------- */
  /*     REFINE ESTIMATE AND TEST */
  /* ----------------------------------------------------------------------- */
  aphi = slatec_zabs_impl_(&phidr, &phidi);
  rs1 += log(aphi);
  if (fabs(rs1) < *elim) {
    goto L100;
  }
L95:
  if (fabs(rs1) > 0.) {
    goto L300;
  }
  /* ----------------------------------------------------------------------- */
  /*     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW */
  /* ----------------------------------------------------------------------- */
  if (*zr < 0.) {
    goto L300;
  }
  *nz = *n;
  i_u_1 = *n;
  for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
    yr[i_u_] = zeror;
    yi[i_u_] = zeroi;
    /* L96: */
  }
  return 0;
/* ----------------------------------------------------------------------- */
/*     FORWARD RECUR FOR REMAINDER OF THE SEQUENCE */
/* ----------------------------------------------------------------------- */
L100:
  s1r = cyr[0];
  s1i = cyi[0];
  s2r = cyr[1];
  s2i = cyi[1];
  c1r = csrr[kflag - 1];
  ascle = bry[kflag - 1];
  i_u_1 = *n;
  for (i_u_ = ib; i_u_ <= i_u_1; ++i_u_) {
    c2r = s2r;
    c2i = s2i;
    s2r = ckr * c2r - cki * c2i + s1r;
    s2i = ckr * c2i + cki * c2r + s1i;
    s1r = c2r;
    s1i = c2i;
    ckr += rzr;
    cki += rzi;
    c2r = s2r * c1r;
    c2i = s2i * c1r;
    yr[i_u_] = c2r;
    yi[i_u_] = c2i;
    if (kflag >= 3) {
      goto L120;
    }
    str = fabs(c2r);
    sti = fabs(c2i);
    c2m = fmax_impl_(str, sti);
    if (c2m <= ascle) {
      goto L120;
    }
    ++kflag;
    ascle = bry[kflag - 1];
    s1r *= c1r;
    s1i *= c1r;
    s2r = c2r;
    s2i = c2i;
    s1r *= cssr[kflag - 1];
    s1i *= cssr[kflag - 1];
    s2r *= cssr[kflag - 1];
    s2i *= cssr[kflag - 1];
    c1r = csrr[kflag - 1];
  L120:;
  }
L160:
  if (*mr == 0) {
    return 0;
  }
  /* ----------------------------------------------------------------------- */
  /*     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0 */
  /* ----------------------------------------------------------------------- */
  *nz = 0;
  fmr = (double)(*mr);
  sgn = -f2c_d_sign_impl_(&pi, &fmr);
  /* ----------------------------------------------------------------------- */
  /*     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP. */
  /* ----------------------------------------------------------------------- */
  csgni = sgn;
  inu = (int)(*fnu);
  fnf = *fnu - inu;
  ifn = inu + *n - 1;
  ang = fnf * sgn;
  cspnr = cos(ang);
  cspni = sin(ang);
  if (ifn % 2 == 0) {
    goto L170;
  }
  cspnr = -cspnr;
  cspni = -cspni;
L170:
  asc = bry[0];
  iuf = 0;
  kk = *n;
  kdflg = 1;
  --ib;
  ic = ib - 1;
  i_u_1 = *n;
  for (k = 1; k <= i_u_1; ++k) {
    fn = *fnu + (kk - 1);
    /* -----------------------------------------------------------------------
     */
    /*     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K */
    /*     FUNCTION ABOVE */
    /* -----------------------------------------------------------------------
     */
    m = 3;
    if (*n > 2) {
      goto L175;
    }
  L172:
    initd = init[j - 1];
    phidr = phir[j - 1];
    phidi = phii[j - 1];
    zet1dr = zeta1r[j - 1];
    zet1di = zeta1i[j - 1];
    zet2dr = zeta2r[j - 1];
    zet2di = zeta2i[j - 1];
    sumdr = sumr[j - 1];
    sumdi = sumi[j - 1];
    m = j;
    j = 3 - j;
    goto L180;
  L175:
    if (kk == *n && ib < *n) {
      goto L180;
    }
    if (kk == ib || kk == ic) {
      goto L172;
    }
    initd = 0;
  L180:
    slatec_zunik_impl_(&zrr, &zri, &fn, &c_u_1, &c_u_0, tol, &initd, &phidr,
                       &phidi, &zet1dr, &zet1di, &zet2dr, &zet2di, &sumdr,
                       &sumdi, &cwrkr[(m << 4) - 16], &cwrki[(m << 4) - 16]);
    if (*kode == 1) {
      goto L200;
    }
    str = zrr + zet2dr;
    sti = zri + zet2di;
    rast = fn / slatec_zabs_impl_(&str, &sti);
    str = str * rast * rast;
    sti = -sti * rast * rast;
    s1r = -zet1dr + str;
    s1i = -zet1di + sti;
    goto L210;
  L200:
    s1r = -zet1dr + zet2dr;
    s1i = -zet1di + zet2di;
  L210:
    /* -----------------------------------------------------------------------
     */
    /*     TEST FOR UNDERFLOW AND OVERFLOW */
    /* -----------------------------------------------------------------------
     */
    rs1 = s1r;
    if (fabs(rs1) > *elim) {
      goto L260;
    }
    if (kdflg == 1) {
      iflag = 2;
    }
    if (fabs(rs1) < *alim) {
      goto L220;
    }
    /* -----------------------------------------------------------------------
     */
    /*     REFINE  TEST AND SCALE */
    /* -----------------------------------------------------------------------
     */
    aphi = slatec_zabs_impl_(&phidr, &phidi);
    rs1 += log(aphi);
    if (fabs(rs1) > *elim) {
      goto L260;
    }
    if (kdflg == 1) {
      iflag = 1;
    }
    if (rs1 < 0.) {
      goto L220;
    }
    if (kdflg == 1) {
      iflag = 3;
    }
  L220:
    str = phidr * sumdr - phidi * sumdi;
    sti = phidr * sumdi + phidi * sumdr;
    s2r = -csgni * sti;
    s2i = csgni * str;
    str = exp(s1r) * cssr[iflag - 1];
    s1r = str * cos(s1i);
    s1i = str * sin(s1i);
    str = s2r * s1r - s2i * s1i;
    s2i = s2r * s1i + s2i * s1r;
    s2r = str;
    if (iflag != 1) {
      goto L230;
    }
    slatec_zuchk_impl_(&s2r, &s2i, &nw, bry, tol);
    if (nw == 0) {
      goto L230;
    }
    s2r = zeror;
    s2i = zeroi;
  L230:
    cyr[kdflg - 1] = s2r;
    cyi[kdflg - 1] = s2i;
    c2r = s2r;
    c2i = s2i;
    s2r *= csrr[iflag - 1];
    s2i *= csrr[iflag - 1];
    /* -----------------------------------------------------------------------
     */
    /*     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N */
    /* -----------------------------------------------------------------------
     */
    s1r = yr[kk];
    s1i = yi[kk];
    if (*kode == 1) {
      goto L250;
    }
    slatec_zs1s2_impl_(&zrr, &zri, &s1r, &s1i, &s2r, &s2i, &nw, &asc, alim,
                       &iuf);
    *nz += nw;
  L250:
    yr[kk] = s1r * cspnr - s1i * cspni + s2r;
    yi[kk] = cspnr * s1i + cspni * s1r + s2i;
    --kk;
    cspnr = -cspnr;
    cspni = -cspni;
    if (c2r != 0. || c2i != 0.) {
      goto L255;
    }
    kdflg = 1;
    goto L270;
  L255:
    if (kdflg == 2) {
      goto L275;
    }
    kdflg = 2;
    goto L270;
  L260:
    if (rs1 > 0.) {
      goto L300;
    }
    s2r = zeror;
    s2i = zeroi;
    goto L230;
  L270:;
  }
  k = *n;
L275:
  il = *n - k;
  if (il == 0) {
    return 0;
  }
  /* ----------------------------------------------------------------------- */
  /*     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE */
  /*     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP */
  /*     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES. */
  /* ----------------------------------------------------------------------- */
  s1r = cyr[0];
  s1i = cyi[0];
  s2r = cyr[1];
  s2i = cyi[1];
  csr = csrr[iflag - 1];
  ascle = bry[iflag - 1];
  fn = (double)(inu + il);
  i_u_1 = il;
  for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
    c2r = s2r;
    c2i = s2i;
    s2r = s1r + (fn + fnf) * (rzr * c2r - rzi * c2i);
    s2i = s1i + (fn + fnf) * (rzr * c2i + rzi * c2r);
    s1r = c2r;
    s1i = c2i;
    fn += -1.;
    c2r = s2r * csr;
    c2i = s2i * csr;
    ckr = c2r;
    cki = c2i;
    c1r = yr[kk];
    c1i = yi[kk];
    if (*kode == 1) {
      goto L280;
    }
    slatec_zs1s2_impl_(&zrr, &zri, &c1r, &c1i, &c2r, &c2i, &nw, &asc, alim,
                       &iuf);
    *nz += nw;
  L280:
    yr[kk] = c1r * cspnr - c1i * cspni + c2r;
    yi[kk] = c1r * cspni + c1i * cspnr + c2i;
    --kk;
    cspnr = -cspnr;
    cspni = -cspni;
    if (iflag >= 3) {
      goto L290;
    }
    c2r = fabs(ckr);
    c2i = fabs(cki);
    c2m = fmax_impl_(c2r, c2i);
    if (c2m <= ascle) {
      goto L290;
    }
    ++iflag;
    ascle = bry[iflag - 1];
    s1r *= csr;
    s1i *= csr;
    s2r = ckr;
    s2i = cki;
    s1r *= cssr[iflag - 1];
    s1i *= cssr[iflag - 1];
    s2r *= cssr[iflag - 1];
    s2i *= cssr[iflag - 1];
    csr = csrr[iflag - 1];
  L290:;
  }
  return 0;
L300:
  *nz = -1;
  return 0;
} /* slatec_zunk1_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZUNK1_IMPL_H */