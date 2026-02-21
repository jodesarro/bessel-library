/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Slatec's ZUNI1 subroutine, originally translated from Fortran to
  C with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZUNI1_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZUNI1_IMPL_H

#include "fmax_impl_.h"
#include "fortran_d1mach_impl_.h"
#include "min_impl_.h"
#include "slatec_zabs_impl_.h"
#include "slatec_zuchk_impl_.h"
#include "slatec_zunik_impl_.h"
#include "slatec_zuoik_impl_.h"
#include <math.h> /* For math operations and constants */

/* zuni1.f -- translated by f2c (version 20100827) */

/* DECK ZUNI1 */
/* Subroutine */
static inline int slatec_zuni1_impl_(double *zr, double *zi, double *fnu,
                                     int *kode, int *n, double *yr, double *yi,
                                     int *nz, int *nlast, double *fnul,
                                     double *tol, double *elim, double *alim) {
  /* Initialized data */
  static int c_u_0 = 0;
  static int c_u_1 = 1;
  static int c_u_2 = 2;
  static double zeror = 0.;
  static double zeroi = 0.;
  static double coner = 1.;

  /* System generated locals */
  int i_u_1;

  /* Local variables */
  static int i_u_, k, m, nd;
  static double fn;
  static int nn, nw;
  static double c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r, cyi[2];
  static int nuf;
  static double bry[3], cyr[2], sti, rzi, str, rzr, aphi, cscl, phii, crsc;
  static double phir;
  static int init;
  static double csrr[3], cssr[3], rast, sumi, sumr;
  static int iflag;
  static double ascle, cwrki[16];
  static double cwrkr[16];
  static double zeta1i, zeta2i, zeta1r, zeta2r;

  /* ***BEGIN PROLOGUE  ZUNI1 */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CUNI1-A, ZUNI1-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC */
  /*     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3. */

  /*     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC */
  /*     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET. */
  /*     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER */
  /*     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL. */
  /*     Y(I)=CZERO FOR I=NLAST+1,N */

  /* ***SEE ALSO  ZBESI, ZBESK */
  /* ***ROUTINES CALLED  D1MACH, ZABS, ZUCHK, ZUNIK, ZUOIK */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZUNI1 */
  /*     COMPLEX CFN,CONE,CRSC,CSCL,CSR,CSS,CWRK,CZERO,C1,C2,PHI,RZ,SUM,S1, */
  /*    *S2,Y,Z,ZETA1,ZETA2 */
  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZUNI1 */
  *nz = 0;
  nd = *n;
  *nlast = 0;
  /* ----------------------------------------------------------------------- */
  /*     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG- */
  /*     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE, */
  /*     EXP(ALIM)=EXP(ELIM)*TOL */
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
  /* ----------------------------------------------------------------------- */
  /*     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER */
  /* ----------------------------------------------------------------------- */
  fn = fmax_impl_(*fnu, 1.);
  init = 0;
  slatec_zunik_impl_(zr, zi, &fn, &c_u_1, &c_u_1, tol, &init, &phir, &phii,
                     &zeta1r, &zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr,
                     cwrki);
  if (*kode == 1) {
    goto L10;
  }
  str = *zr + zeta2r;
  sti = *zi + zeta2i;
  rast = fn / slatec_zabs_impl_(&str, &sti);
  str = str * rast * rast;
  sti = -sti * rast * rast;
  s1r = -zeta1r + str;
  s1i = -zeta1i + sti;
  goto L20;
L10:
  s1r = -zeta1r + zeta2r;
  s1i = -zeta1i + zeta2i;
L20:
  rs1 = s1r;
  if (fabs(rs1) > *elim) {
    goto L130;
  }
L30:
  nn = min_impl_(2, nd);
  i_u_1 = nn;
  for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
    fn = *fnu + (nd - i_u_);
    init = 0;
    slatec_zunik_impl_(zr, zi, &fn, &c_u_1, &c_u_0, tol, &init, &phir, &phii,
                       &zeta1r, &zeta1i, &zeta2r, &zeta2i, &sumr, &sumi, cwrkr,
                       cwrki);
    if (*kode == 1) {
      goto L40;
    }
    str = *zr + zeta2r;
    sti = *zi + zeta2i;
    rast = fn / slatec_zabs_impl_(&str, &sti);
    str = str * rast * rast;
    sti = -sti * rast * rast;
    s1r = -zeta1r + str;
    s1i = -zeta1i + sti + *zi;
    goto L50;
  L40:
    s1r = -zeta1r + zeta2r;
    s1i = -zeta1i + zeta2i;
  L50:
    /* -----------------------------------------------------------------------
     */
    /*     TEST FOR UNDERFLOW AND OVERFLOW */
    /* -----------------------------------------------------------------------
     */
    rs1 = s1r;
    if (fabs(rs1) > *elim) {
      goto L110;
    }
    if (i_u_ == 1) {
      iflag = 2;
    }
    if (fabs(rs1) < *alim) {
      goto L60;
    }
    /* -----------------------------------------------------------------------
     */
    /*     REFINE  TEST AND SCALE */
    /* -----------------------------------------------------------------------
     */
    aphi = slatec_zabs_impl_(&phir, &phii);
    rs1 += log(aphi);
    if (fabs(rs1) > *elim) {
      goto L110;
    }
    if (i_u_ == 1) {
      iflag = 1;
    }
    if (rs1 < 0.) {
      goto L60;
    }
    if (i_u_ == 1) {
      iflag = 3;
    }
  L60:
    /* -----------------------------------------------------------------------
     */
    /*     SCALE S1 IF ABS(S1).LT.ASCLE */
    /* -----------------------------------------------------------------------
     */
    s2r = phir * sumr - phii * sumi;
    s2i = phir * sumi + phii * sumr;
    str = exp(s1r) * cssr[iflag - 1];
    s1r = str * cos(s1i);
    s1i = str * sin(s1i);
    str = s2r * s1r - s2i * s1i;
    s2i = s2r * s1i + s2i * s1r;
    s2r = str;
    if (iflag != 1) {
      goto L70;
    }
    slatec_zuchk_impl_(&s2r, &s2i, &nw, bry, tol);
    if (nw != 0) {
      goto L110;
    }
  L70:
    cyr[i_u_ - 1] = s2r;
    cyi[i_u_ - 1] = s2i;
    m = nd - i_u_ + 1;
    yr[m] = s2r * csrr[iflag - 1];
    yi[m] = s2i * csrr[iflag - 1];
    /* L80: */
  }
  if (nd <= 2) {
    goto L100;
  }
  rast = 1. / slatec_zabs_impl_(zr, zi);
  str = *zr * rast;
  sti = -(*zi) * rast;
  rzr = (str + str) * rast;
  rzi = (sti + sti) * rast;
  bry[1] = 1. / bry[0];
  bry[2] = fortran_d1mach_impl_(&c_u_2);
  s1r = cyr[0];
  s1i = cyi[0];
  s2r = cyr[1];
  s2i = cyi[1];
  c1r = csrr[iflag - 1];
  ascle = bry[iflag - 1];
  k = nd - 2;
  fn = (double)k;
  i_u_1 = nd;
  for (i_u_ = 3; i_u_ <= i_u_1; ++i_u_) {
    c2r = s2r;
    c2i = s2i;
    s2r = s1r + (*fnu + fn) * (rzr * c2r - rzi * c2i);
    s2i = s1i + (*fnu + fn) * (rzr * c2i + rzi * c2r);
    s1r = c2r;
    s1i = c2i;
    c2r = s2r * c1r;
    c2i = s2i * c1r;
    yr[k] = c2r;
    yi[k] = c2i;
    --k;
    fn += -1.;
    if (iflag >= 3) {
      goto L90;
    }
    str = fabs(c2r);
    sti = fabs(c2i);
    c2m = fmax_impl_(str, sti);
    if (c2m <= ascle) {
      goto L90;
    }
    ++iflag;
    ascle = bry[iflag - 1];
    s1r *= c1r;
    s1i *= c1r;
    s2r = c2r;
    s2i = c2i;
    s1r *= cssr[iflag - 1];
    s1i *= cssr[iflag - 1];
    s2r *= cssr[iflag - 1];
    s2i *= cssr[iflag - 1];
    c1r = csrr[iflag - 1];
  L90:;
  }
L100:
  return 0;
/* ----------------------------------------------------------------------- */
/*     SET UNDERFLOW AND UPDATE PARAMETERS */
/* ----------------------------------------------------------------------- */
L110:
  if (rs1 > 0.) {
    goto L120;
  }
  yr[nd] = zeror;
  yi[nd] = zeroi;
  ++(*nz);
  --nd;
  if (nd == 0) {
    goto L100;
  }
  slatec_zuoik_impl_(zr, zi, fnu, kode, &c_u_1, &nd, &yr[1], &yi[1], &nuf, tol,
                     elim, alim);
  if (nuf < 0) {
    goto L120;
  }
  nd -= nuf;
  *nz += nuf;
  if (nd == 0) {
    goto L100;
  }
  fn = *fnu + (nd - 1);
  if (fn >= *fnul) {
    goto L30;
  }
  *nlast = nd;
  return 0;
L120:
  *nz = -1;
  return 0;
L130:
  if (rs1 > 0.) {
    goto L120;
  }
  *nz = *n;
  i_u_1 = *n;
  for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
    yr[i_u_] = zeror;
    yi[i_u_] = zeroi;
    /* L140: */
  }
  return 0;
} /* slatec_zuni1_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZUNI1_IMPL_H */