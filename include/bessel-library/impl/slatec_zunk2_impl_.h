/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Slatec's ZUNK2 subroutine, originally translated from Fortran to
  C with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZUNK2_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZUNK2_IMPL_H

#include "f2c_d_sign_impl_.h"
#include "fmax_impl_.h"
#include "fortran_d1mach_impl_.h"
#include "slatec_zabs_impl_.h"
#include "slatec_zairy_impl_.h"
#include "slatec_zs1s2_impl_.h"
#include "slatec_zuchk_impl_.h"
#include "slatec_zunhj_impl_.h"
#include <math.h> /* For math operations and constants */

/* zunk2.f -- translated by f2c (version 20100827) */

/* DECK ZUNK2 */
/* Subroutine */
static inline int slatec_zunk2_impl_(double *zr, double *zi, double *fnu,
                                     int *kode, int *mr, int *n, double *yr,
                                     double *yi, int *nz, double *tol,
                                     double *elim, double *alim) {
  /* Initialized data */
  int c_u_0 = 0;
  int c_u_1 = 1;
  int c_u_2 = 2;
  double zeror = 0.;
  double aic = 1.26551212348464539;
  double cipr[4] = {1., 0., -1., 0.};
  double cipi[4] = {0., -1., 0., 1.};
  double zeroi = 0.;
  double coner = 1.;
  double cr1r = 1.;
  double cr1i = 1.73205080756887729;
  double cr2r = -.5;
  double cr2i = -.866025403784438647;
  double hpi = 1.57079632679489662;
  double pi = 3.14159265358979324;

  /* System generated locals */
  int i_u_1;

  /* Local variables */
  int i_u_, j, k, ib, ic;
  double fn;
  int il, kk, in, nw;
  double yy, c1i, c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r, aii, ang, asc,
      car, cki, fnf;
  int nai;
  double air;
  int ifn;
  double csi, ckr;
  int iuf;
  double cyi[2], fmr, sar, csr, sgn, zbi;
  int inu;
  double bry[3], cyr[2], pti, sti, zbr, zni, rzi, ptr, zri, str, znr, rzr, zrr,
      daii, aarg;
  int ndai;
  double dair, aphi, argi[2], cscl, phii[2], crsc, argr[2];
  int idum;
  double phir[2], csrr[3], cssr[3], rast, razr;
  int iflag, kflag;
  double argdi, ascle;
  int kdflg;
  double phidi, argdr;
  int ipard;
  double csgni, phidr, cspni, asumi[2], bsumi[2];
  double cspnr, asumr[2], bsumr[2];
  double zeta1i[2], zeta2i[2], zet1di, zet2di, zeta1r[2], zeta2r[2], zet1dr,
      zet2dr, asumdi, bsumdi, asumdr, bsumdr;

  /* ***BEGIN PROLOGUE  ZUNK2 */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CUNK2-A, ZUNK2-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE */
  /*     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE */
  /*     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN) */
  /*     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR */
  /*     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT */
  /*     HALF PLANE OR ZR=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC- */
  /*     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION. */
  /*     NZ=-1 MEANS AN OVERFLOW WILL OCCUR */

  /* ***SEE ALSO  ZBESK */
  /* ***ROUTINES CALLED  D1MACH, ZABS, ZAIRY, ZS1S2, ZUCHK, ZUNHJ */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZUNK2 */
  /*     COMPLEX AI,ARG,ARGD,ASUM,ASUMD,BSUM,BSUMD,CFN,CI,CIP,CK,CONE,CRSC, */
  /*    *CR1,CR2,CS,CSCL,CSGN,CSPN,CSR,CSS,CY,CZERO,C1,C2,DAI,PHI,PHID,RZ, */
  /*    *S1,S2,Y,Z,ZB,ZETA1,ZETA1D,ZETA2,ZETA2D,ZN,ZR */
  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZUNK2 */
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
  yy = zri;
  znr = zri;
  zni = -zrr;
  zbr = zrr;
  zbi = zri;
  inu = (int)(*fnu);
  fnf = *fnu - inu;
  ang = -hpi * fnf;
  car = cos(ang);
  sar = sin(ang);
  c2r = hpi * sar;
  c2i = -hpi * car;
  kk = inu % 4 + 1;
  str = c2r * cipr[kk - 1] - c2i * cipi[kk - 1];
  sti = c2r * cipi[kk - 1] + c2i * cipr[kk - 1];
  csr = cr1r * str - cr1i * sti;
  csi = cr1r * sti + cr1i * str;
  if (yy > 0.) {
    goto L20;
  }
  znr = -znr;
  zbi = -zbi;
L20:
  /* ----------------------------------------------------------------------- */
  /*     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST */
  /*     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY */
  /*     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS */
  /* ----------------------------------------------------------------------- */
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
    slatec_zunhj_impl_(&znr, &zni, &fn, &c_u_0, tol, &phir[j - 1], &phii[j - 1],
                       &argr[j - 1], &argi[j - 1], &zeta1r[j - 1],
                       &zeta1i[j - 1], &zeta2r[j - 1], &zeta2i[j - 1],
                       &asumr[j - 1], &asumi[j - 1], &bsumr[j - 1],
                       &bsumi[j - 1]);
    if (*kode == 1) {
      goto L30;
    }
    str = zbr + zeta2r[j - 1];
    sti = zbi + zeta2i[j - 1];
    rast = fn / slatec_zabs_impl_(&str, &sti);
    str = str * rast * rast;
    sti = -sti * rast * rast;
    s1r = zeta1r[j - 1] - str;
    s1i = zeta1i[j - 1] - sti;
    goto L40;
  L30:
    s1r = zeta1r[j - 1] - zeta2r[j - 1];
    s1i = zeta1i[j - 1] - zeta2i[j - 1];
  L40:
    /* -----------------------------------------------------------------------
     */
    /*     TEST FOR UNDERFLOW AND OVERFLOW */
    /* -----------------------------------------------------------------------
     */
    rs1 = s1r;
    if (fabs(rs1) > *elim) {
      goto L70;
    }
    if (kdflg == 1) {
      kflag = 2;
    }
    if (fabs(rs1) < *alim) {
      goto L50;
    }
    /* -----------------------------------------------------------------------
     */
    /*     REFINE  TEST AND SCALE */
    /* -----------------------------------------------------------------------
     */
    aphi = slatec_zabs_impl_(&phir[j - 1], &phii[j - 1]);
    aarg = slatec_zabs_impl_(&argr[j - 1], &argi[j - 1]);
    rs1 = rs1 + log(aphi) - log(aarg) * .25 - aic;
    if (fabs(rs1) > *elim) {
      goto L70;
    }
    if (kdflg == 1) {
      kflag = 1;
    }
    if (rs1 < 0.) {
      goto L50;
    }
    if (kdflg == 1) {
      kflag = 3;
    }
  L50:
    /* -----------------------------------------------------------------------
     */
    /*     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR */
    /*     EXPONENT EXTREMES */
    /* -----------------------------------------------------------------------
     */
    c2r = argr[j - 1] * cr2r - argi[j - 1] * cr2i;
    c2i = argr[j - 1] * cr2i + argi[j - 1] * cr2r;
    slatec_zairy_impl_(&c2r, &c2i, &c_u_0, &c_u_2, &air, &aii, &nai, &idum);
    slatec_zairy_impl_(&c2r, &c2i, &c_u_1, &c_u_2, &dair, &daii, &ndai, &idum);
    str = dair * bsumr[j - 1] - daii * bsumi[j - 1];
    sti = dair * bsumi[j - 1] + daii * bsumr[j - 1];
    ptr = str * cr2r - sti * cr2i;
    pti = str * cr2i + sti * cr2r;
    str = ptr + (air * asumr[j - 1] - aii * asumi[j - 1]);
    sti = pti + (air * asumi[j - 1] + aii * asumr[j - 1]);
    ptr = str * phir[j - 1] - sti * phii[j - 1];
    pti = str * phii[j - 1] + sti * phir[j - 1];
    s2r = ptr * csr - pti * csi;
    s2i = ptr * csi + pti * csr;
    str = exp(s1r) * cssr[kflag - 1];
    s1r = str * cos(s1i);
    s1i = str * sin(s1i);
    str = s2r * s1r - s2i * s1i;
    s2i = s1r * s2i + s2r * s1i;
    s2r = str;
    if (kflag != 1) {
      goto L60;
    }
    slatec_zuchk_impl_(&s2r, &s2i, &nw, bry, tol);
    if (nw != 0) {
      goto L70;
    }
  L60:
    if (yy <= 0.) {
      s2i = -s2i;
    }
    cyr[kdflg - 1] = s2r;
    cyi[kdflg - 1] = s2i;
    yr[i_u_] = s2r * csrr[kflag - 1];
    yi[i_u_] = s2i * csrr[kflag - 1];
    str = csi;
    csi = -csr;
    csr = str;
    if (kdflg == 2) {
      goto L85;
    }
    kdflg = 2;
    goto L80;
  L70:
    if (rs1 > 0.) {
      goto L320;
    }
    /* -----------------------------------------------------------------------
     */
    /*     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW */
    /* -----------------------------------------------------------------------
     */
    if (*zr < 0.) {
      goto L320;
    }
    kdflg = 1;
    yr[i_u_] = zeror;
    yi[i_u_] = zeroi;
    ++(*nz);
    str = csi;
    csi = -csr;
    csr = str;
    if (i_u_ == 1) {
      goto L80;
    }
    if (yr[i_u_ - 1] == zeror && yi[i_u_ - 1] == zeroi) {
      goto L80;
    }
    yr[i_u_ - 1] = zeror;
    yi[i_u_ - 1] = zeroi;
    ++(*nz);
  L80:;
  }
  i_u_ = *n;
L85:
  razr = 1. / slatec_zabs_impl_(&zrr, &zri);
  str = zrr * razr;
  sti = -zri * razr;
  rzr = (str + str) * razr;
  rzi = (sti + sti) * razr;
  ckr = fn * rzr;
  cki = fn * rzi;
  ib = i_u_ + 1;
  if (*n < ib) {
    goto L180;
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
  slatec_zunhj_impl_(&znr, &zni, &fn, &ipard, tol, &phidr, &phidi, &argdr,
                     &argdi, &zet1dr, &zet1di, &zet2dr, &zet2di, &asumdr,
                     &asumdi, &bsumdr, &bsumdi);
  if (*kode == 1) {
    goto L90;
  }
  str = zbr + zet2dr;
  sti = zbi + zet2di;
  rast = fn / slatec_zabs_impl_(&str, &sti);
  str = str * rast * rast;
  sti = -sti * rast * rast;
  s1r = zet1dr - str;
  s1i = zet1di - sti;
  goto L100;
L90:
  s1r = zet1dr - zet2dr;
  s1i = zet1di - zet2di;
L100:
  rs1 = s1r;
  if (fabs(rs1) > *elim) {
    goto L105;
  }
  if (fabs(rs1) < *alim) {
    goto L120;
  }
  /* ----------------------------------------------------------------------- */
  /*     REFINE ESTIMATE AND TEST */
  /* ----------------------------------------------------------------------- */
  aphi = slatec_zabs_impl_(&phidr, &phidi);
  rs1 += log(aphi);
  if (fabs(rs1) < *elim) {
    goto L120;
  }
L105:
  if (rs1 > 0.) {
    goto L320;
  }
  /* ----------------------------------------------------------------------- */
  /*     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW */
  /* ----------------------------------------------------------------------- */
  if (*zr < 0.) {
    goto L320;
  }
  *nz = *n;
  i_u_1 = *n;
  for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
    yr[i_u_] = zeror;
    yi[i_u_] = zeroi;
    /* L106: */
  }
  return 0;
L120:
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
      goto L130;
    }
    str = fabs(c2r);
    sti = fabs(c2i);
    c2m = fmax_impl_(str, sti);
    if (c2m <= ascle) {
      goto L130;
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
  L130:;
  }
L180:
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
  if (yy <= 0.) {
    csgni = -csgni;
  }
  ifn = inu + *n - 1;
  ang = fnf * sgn;
  cspnr = cos(ang);
  cspni = sin(ang);
  if (ifn % 2 == 0) {
    goto L190;
  }
  cspnr = -cspnr;
  cspni = -cspni;
L190:
  /* ----------------------------------------------------------------------- */
  /*     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS */
  /*     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST */
  /*     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY */
  /*     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS */
  /* ----------------------------------------------------------------------- */
  csr = sar * csgni;
  csi = car * csgni;
  in = ifn % 4 + 1;
  c2r = cipr[in - 1];
  c2i = cipi[in - 1];
  str = csr * c2r + csi * c2i;
  csi = -csr * c2i + csi * c2r;
  csr = str;
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
    if (*n > 2) {
      goto L175;
    }
  L172:
    phidr = phir[j - 1];
    phidi = phii[j - 1];
    argdr = argr[j - 1];
    argdi = argi[j - 1];
    zet1dr = zeta1r[j - 1];
    zet1di = zeta1i[j - 1];
    zet2dr = zeta2r[j - 1];
    zet2di = zeta2i[j - 1];
    asumdr = asumr[j - 1];
    asumdi = asumi[j - 1];
    bsumdr = bsumr[j - 1];
    bsumdi = bsumi[j - 1];
    j = 3 - j;
    goto L210;
  L175:
    if (kk == *n && ib < *n) {
      goto L210;
    }
    if (kk == ib || kk == ic) {
      goto L172;
    }
    slatec_zunhj_impl_(&znr, &zni, &fn, &c_u_0, tol, &phidr, &phidi, &argdr,
                       &argdi, &zet1dr, &zet1di, &zet2dr, &zet2di, &asumdr,
                       &asumdi, &bsumdr, &bsumdi);
  L210:
    if (*kode == 1) {
      goto L220;
    }
    str = zbr + zet2dr;
    sti = zbi + zet2di;
    rast = fn / slatec_zabs_impl_(&str, &sti);
    str = str * rast * rast;
    sti = -sti * rast * rast;
    s1r = -zet1dr + str;
    s1i = -zet1di + sti;
    goto L230;
  L220:
    s1r = -zet1dr + zet2dr;
    s1i = -zet1di + zet2di;
  L230:
    /* -----------------------------------------------------------------------
     */
    /*     TEST FOR UNDERFLOW AND OVERFLOW */
    /* -----------------------------------------------------------------------
     */
    rs1 = s1r;
    if (fabs(rs1) > *elim) {
      goto L280;
    }
    if (kdflg == 1) {
      iflag = 2;
    }
    if (fabs(rs1) < *alim) {
      goto L240;
    }
    /* -----------------------------------------------------------------------
     */
    /*     REFINE  TEST AND SCALE */
    /* -----------------------------------------------------------------------
     */
    aphi = slatec_zabs_impl_(&phidr, &phidi);
    aarg = slatec_zabs_impl_(&argdr, &argdi);
    rs1 = rs1 + log(aphi) - log(aarg) * .25 - aic;
    if (fabs(rs1) > *elim) {
      goto L280;
    }
    if (kdflg == 1) {
      iflag = 1;
    }
    if (rs1 < 0.) {
      goto L240;
    }
    if (kdflg == 1) {
      iflag = 3;
    }
  L240:
    slatec_zairy_impl_(&argdr, &argdi, &c_u_0, &c_u_2, &air, &aii, &nai, &idum);
    slatec_zairy_impl_(&argdr, &argdi, &c_u_1, &c_u_2, &dair, &daii, &ndai,
                       &idum);
    str = dair * bsumdr - daii * bsumdi;
    sti = dair * bsumdi + daii * bsumdr;
    str += air * asumdr - aii * asumdi;
    sti += air * asumdi + aii * asumdr;
    ptr = str * phidr - sti * phidi;
    pti = str * phidi + sti * phidr;
    s2r = ptr * csr - pti * csi;
    s2i = ptr * csi + pti * csr;
    str = exp(s1r) * cssr[iflag - 1];
    s1r = str * cos(s1i);
    s1i = str * sin(s1i);
    str = s2r * s1r - s2i * s1i;
    s2i = s2r * s1i + s2i * s1r;
    s2r = str;
    if (iflag != 1) {
      goto L250;
    }
    slatec_zuchk_impl_(&s2r, &s2i, &nw, bry, tol);
    if (nw == 0) {
      goto L250;
    }
    s2r = zeror;
    s2i = zeroi;
  L250:
    if (yy <= 0.) {
      s2i = -s2i;
    }
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
      goto L270;
    }
    slatec_zs1s2_impl_(&zrr, &zri, &s1r, &s1i, &s2r, &s2i, &nw, &asc, alim,
                       &iuf);
    *nz += nw;
  L270:
    yr[kk] = s1r * cspnr - s1i * cspni + s2r;
    yi[kk] = s1r * cspni + s1i * cspnr + s2i;
    --kk;
    cspnr = -cspnr;
    cspni = -cspni;
    str = csi;
    csi = -csr;
    csr = str;
    if (c2r != 0. || c2i != 0.) {
      goto L255;
    }
    kdflg = 1;
    goto L290;
  L255:
    if (kdflg == 2) {
      goto L295;
    }
    kdflg = 2;
    goto L290;
  L280:
    if (rs1 > 0.) {
      goto L320;
    }
    s2r = zeror;
    s2i = zeroi;
    goto L250;
  L290:;
  }
  k = *n;
L295:
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
      goto L300;
    }
    slatec_zs1s2_impl_(&zrr, &zri, &c1r, &c1i, &c2r, &c2i, &nw, &asc, alim,
                       &iuf);
    *nz += nw;
  L300:
    yr[kk] = c1r * cspnr - c1i * cspni + c2r;
    yi[kk] = c1r * cspni + c1i * cspnr + c2i;
    --kk;
    cspnr = -cspnr;
    cspni = -cspni;
    if (iflag >= 3) {
      goto L310;
    }
    c2r = fabs(ckr);
    c2i = fabs(cki);
    c2m = fmax_impl_(c2r, c2i);
    if (c2m <= ascle) {
      goto L310;
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
  L310:;
  }
  return 0;
L320:
  *nz = -1;
  return 0;
} /* slatec_zunk2_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZUNK2_IMPL_H */