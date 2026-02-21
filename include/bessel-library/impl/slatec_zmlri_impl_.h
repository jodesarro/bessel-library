/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Slatec's ZMLRI subroutine, originally translated from Fortran to
  C with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZMLRI_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZMLRI_IMPL_H

#include "fmin_impl_.h"
#include "fortran_d1mach_impl_.h"
#include "max_impl_.h"
#include "slatec_dgamln_impl_.h"
#include "slatec_zabs_impl_.h"
#include "slatec_zexp_impl_.h"
#include "slatec_zlog_impl_.h"
#include "slatec_zmlt_impl_.h"
#include <math.h> /* For math operations and constants */

/* zmlri.f -- translated by f2c (version 20100827) */

/* DECK ZMLRI */
/* Subroutine */
static inline int slatec_zmlri_impl_(double *zr, double *zi, double *fnu,
                                     int *kode, int *n, double *yr, double *yi,
                                     int *nz, double *tol) {
  /* Initialized data */
  static int c_u_1 = 1;
  static double zeror = 0.;
  static double zeroi = 0.;
  static double coner = 1.;
  static double conei = 0.;

  /* System generated locals */
  int i_u_1, i_u_2;
  double d_u_1, d_u_2, d_u_3;

  /* Local variables */
  static int i_u_, k, m;
  static double ak, bk, ap, at;
  static int kk, km;
  static double az, p1i, p2i, p1r, p2r, ack, cki, fnf, fkk, ckr;
  static int iaz;
  static double rho;
  static int inu;
  static double pti, raz, sti, rzi, ptr, str, tst, rzr, rho2, flam, fkap, scle,
      tfnf;
  static int idum;
  static int ifnu;
  static double sumi, sumr;
  static int itime;
  static double cnormi, cnormr;

  /* ***BEGIN PROLOGUE  ZMLRI */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CMLRI-A, ZMLRI-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE */
  /*     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES. */

  /* ***SEE ALSO  ZBESI, ZBESK */
  /* ***ROUTINES CALLED  D1MACH, DGAMLN, ZABS, ZEXP, ZLOG, ZMLT */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /*   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC) */
  /* ***END PROLOGUE  ZMLRI */
  /*     COMPLEX CK,CNORM,CONE,CTWO,CZERO,PT,P1,P2,RZ,SUM,Y,Z */
  /* Parameter adjustments */
  --yi;
  --yr;

  /* Function Body */
  /* ***FIRST EXECUTABLE STATEMENT  ZMLRI */
  scle = fortran_d1mach_impl_(&c_u_1) / *tol;
  *nz = 0;
  az = slatec_zabs_impl_(zr, zi);
  iaz = (int)az;
  ifnu = (int)(*fnu);
  inu = ifnu + *n - 1;
  at = iaz + 1.;
  raz = 1. / az;
  str = *zr * raz;
  sti = -(*zi) * raz;
  ckr = str * at * raz;
  cki = sti * at * raz;
  rzr = (str + str) * raz;
  rzi = (sti + sti) * raz;
  p1r = zeror;
  p1i = zeroi;
  p2r = coner;
  p2i = conei;
  ack = (at + 1.) * raz;
  rho = ack + sqrt(ack * ack - 1.);
  rho2 = rho * rho;
  tst = (rho2 + rho2) / ((rho2 - 1.) * (rho - 1.));
  tst /= *tol;
  /* ----------------------------------------------------------------------- */
  /*     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES */
  /* ----------------------------------------------------------------------- */
  ak = at;
  for (i_u_ = 1; i_u_ <= 80; ++i_u_) {
    ptr = p2r;
    pti = p2i;
    p2r = p1r - (ckr * ptr - cki * pti);
    p2i = p1i - (cki * ptr + ckr * pti);
    p1r = ptr;
    p1i = pti;
    ckr += rzr;
    cki += rzi;
    ap = slatec_zabs_impl_(&p2r, &p2i);
    if (ap > tst * ak * ak) {
      goto L20;
    }
    ak += 1.;
    /* L10: */
  }
  goto L110;
L20:
  ++i_u_;
  k = 0;
  if (inu < iaz) {
    goto L40;
  }
  /* ----------------------------------------------------------------------- */
  /*     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS */
  /* ----------------------------------------------------------------------- */
  p1r = zeror;
  p1i = zeroi;
  p2r = coner;
  p2i = conei;
  at = inu + 1.;
  str = *zr * raz;
  sti = -(*zi) * raz;
  ckr = str * at * raz;
  cki = sti * at * raz;
  ack = at * raz;
  tst = sqrt(ack / *tol);
  itime = 1;
  for (k = 1; k <= 80; ++k) {
    ptr = p2r;
    pti = p2i;
    p2r = p1r - (ckr * ptr - cki * pti);
    p2i = p1i - (ckr * pti + cki * ptr);
    p1r = ptr;
    p1i = pti;
    ckr += rzr;
    cki += rzi;
    ap = slatec_zabs_impl_(&p2r, &p2i);
    if (ap < tst) {
      goto L30;
    }
    if (itime == 2) {
      goto L40;
    }
    ack = slatec_zabs_impl_(&ckr, &cki);
    flam = ack + sqrt(ack * ack - 1.);
    fkap = ap / slatec_zabs_impl_(&p1r, &p1i);
    rho = fmin_impl_(flam, fkap);
    tst *= sqrt(rho / (rho * rho - 1.));
    itime = 2;
  L30:;
  }
  goto L110;
L40:
  /* ----------------------------------------------------------------------- */
  /*     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION */
  /* ----------------------------------------------------------------------- */
  ++k;
  /* Computing MAX */
  i_u_1 = i_u_ + iaz, i_u_2 = k + inu;
  kk = max_impl_(i_u_1, i_u_2);
  fkk = (double)kk;
  p1r = zeror;
  p1i = zeroi;
  /* ----------------------------------------------------------------------- */
  /*     SCALE P2 AND SUM BY SCLE */
  /* ----------------------------------------------------------------------- */
  p2r = scle;
  p2i = zeroi;
  fnf = *fnu - ifnu;
  tfnf = fnf + fnf;
  d_u_1 = fkk + tfnf + 1.;
  d_u_2 = fkk + 1.;
  d_u_3 = tfnf + 1.;
  bk = slatec_dgamln_impl_(&d_u_1, &idum) - slatec_dgamln_impl_(&d_u_2, &idum) -
       slatec_dgamln_impl_(&d_u_3, &idum);
  bk = exp(bk);
  sumr = zeror;
  sumi = zeroi;
  km = kk - inu;
  i_u_1 = km;
  for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
    ptr = p2r;
    pti = p2i;
    p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
    p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
    p1r = ptr;
    p1i = pti;
    ak = 1. - tfnf / (fkk + tfnf);
    ack = bk * ak;
    sumr += (ack + bk) * p1r;
    sumi += (ack + bk) * p1i;
    bk = ack;
    fkk += -1.;
    /* L50: */
  }
  yr[*n] = p2r;
  yi[*n] = p2i;
  if (*n == 1) {
    goto L70;
  }
  i_u_1 = *n;
  for (i_u_ = 2; i_u_ <= i_u_1; ++i_u_) {
    ptr = p2r;
    pti = p2i;
    p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
    p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
    p1r = ptr;
    p1i = pti;
    ak = 1. - tfnf / (fkk + tfnf);
    ack = bk * ak;
    sumr += (ack + bk) * p1r;
    sumi += (ack + bk) * p1i;
    bk = ack;
    fkk += -1.;
    m = *n - i_u_ + 1;
    yr[m] = p2r;
    yi[m] = p2i;
    /* L60: */
  }
L70:
  if (ifnu <= 0) {
    goto L90;
  }
  i_u_1 = ifnu;
  for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
    ptr = p2r;
    pti = p2i;
    p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
    p2i = p1i + (fkk + fnf) * (rzr * pti + rzi * ptr);
    p1r = ptr;
    p1i = pti;
    ak = 1. - tfnf / (fkk + tfnf);
    ack = bk * ak;
    sumr += (ack + bk) * p1r;
    sumi += (ack + bk) * p1i;
    bk = ack;
    fkk += -1.;
    /* L80: */
  }
L90:
  ptr = *zr;
  pti = *zi;
  if (*kode == 2) {
    ptr = zeror;
  }
  slatec_zlog_impl_(&rzr, &rzi, &str, &sti, &idum);
  p1r = -fnf * str + ptr;
  p1i = -fnf * sti + pti;
  d_u_1 = fnf + 1.;
  ap = slatec_dgamln_impl_(&d_u_1, &idum);
  ptr = p1r - ap;
  pti = p1i;
  /* ----------------------------------------------------------------------- */
  /*     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW */
  /*     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES */
  /* ----------------------------------------------------------------------- */
  p2r += sumr;
  p2i += sumi;
  ap = slatec_zabs_impl_(&p2r, &p2i);
  p1r = 1. / ap;
  slatec_zexp_impl_(&ptr, &pti, &str, &sti);
  ckr = str * p1r;
  cki = sti * p1r;
  ptr = p2r * p1r;
  pti = -p2i * p1r;
  slatec_zmlt_impl_(&ckr, &cki, &ptr, &pti, &cnormr, &cnormi);
  i_u_1 = *n;
  for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
    str = yr[i_u_] * cnormr - yi[i_u_] * cnormi;
    yi[i_u_] = yr[i_u_] * cnormi + yi[i_u_] * cnormr;
    yr[i_u_] = str;
    /* L100: */
  }
  return 0;
L110:
  *nz = -2;
  return 0;
} /* slatec_zmlri_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZMLRI_IMPL_H */