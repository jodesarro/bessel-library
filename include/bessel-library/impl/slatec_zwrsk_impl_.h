/*
  Bessel Library: A C library with routines for computing Bessel functions

  File: include/bessel-library/impl/slatec_zwrsk_impl_.h
  Language standards: C99
  References: include/bessel-library/references.txt
  License: include/bessel-library/license.txt
  Repository: <https://github.com/jodesarro/bessel-library>

  Description: Slatec's ZWRSK subroutine, originally translated from Fortran to
  C with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZWRSK_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZWRSK_IMPL_H

#include "fortran_d1mach_impl_.h"
#include "slatec_zabs_impl_.h"
#include "slatec_zbknu_impl_.h"
#include "slatec_zrati_impl_.h"
#include <math.h> /* For math operations and constants */

/* zwrsk.f -- translated by f2c (version 20100827) */

/* DECK ZWRSK */
/* Subroutine */
static inline int slatec_zwrsk_impl_(double *zrr, double *zri, double *fnu,
                                     int *kode, int *n, double *yr, double *yi,
                                     int *nz, double *cwr, double *cwi,
                                     double *tol, double *elim, double *alim) {
  /* Initialized data */
  static int c_u_1 = 1;
  static int c_u_2 = 2;

  /* System generated locals */
  int i_u_1;

  /* Local variables */
  static int i_u_, nw;
  static double c1i, c2i, c1r, c2r, act, acw, cti, ctr, pti, sti, ptr, str,
      ract;
  static double ascle, csclr, cinui, cinur;

  /* ***BEGIN PROLOGUE  ZWRSK */
  /* ***SUBSIDIARY */
  /* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
  /* ***LIBRARY   SLATEC */
  /* ***TYPE      ALL (CWRSK-A, ZWRSK-A) */
  /* ***AUTHOR  Amos, D. E., (SNL) */
  /* ***DESCRIPTION */

  /*     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY */
  /*     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN */

  /* ***SEE ALSO  ZBESI, ZBESK */
  /* ***ROUTINES CALLED  D1MACH, ZABS, ZBKNU, ZRATI */
  /* ***REVISION HISTORY  (YYMMDD) */
  /*   830501  DATE WRITTEN */
  /*   910415  Prologue converted to Version 4.0 format.  (BAB) */
  /* ***END PROLOGUE  ZWRSK */
  /*     COMPLEX CINU,CSCL,CT,CW,C1,C2,RCT,ST,Y,ZR */
  /* ***FIRST EXECUTABLE STATEMENT  ZWRSK */
  /* ----------------------------------------------------------------------- */
  /*     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS */
  /*     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE */
  /*     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU. */
  /* ----------------------------------------------------------------------- */

  /* Parameter adjustments */
  --yi;
  --yr;
  --cwr;
  --cwi;

  /* Function Body */
  *nz = 0;
  slatec_zbknu_impl_(zrr, zri, fnu, kode, &c_u_2, &cwr[1], &cwi[1], &nw, tol,
                     elim, alim);
  if (nw != 0) {
    goto L50;
  }
  slatec_zrati_impl_(zrr, zri, fnu, n, &yr[1], &yi[1], tol);
  /* ----------------------------------------------------------------------- */
  /*     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z), */
  /*     R(FNU+J-1,Z)=Y(J),  J=1,...,N */
  /* ----------------------------------------------------------------------- */
  cinur = 1.;
  cinui = 0.;
  if (*kode == 1) {
    goto L10;
  }
  cinur = cos(*zri);
  cinui = sin(*zri);
L10:
  /* ----------------------------------------------------------------------- */
  /*     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH */
  /*     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE */
  /*     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT */
  /*     THE RESULT IS ON SCALE. */
  /* ----------------------------------------------------------------------- */
  acw = slatec_zabs_impl_(&cwr[2], &cwi[2]);
  ascle = fortran_d1mach_impl_(&c_u_1) * 1e3 / *tol;
  csclr = 1.;
  if (acw > ascle) {
    goto L20;
  }
  csclr = 1. / *tol;
  goto L30;
L20:
  ascle = 1. / ascle;
  if (acw < ascle) {
    goto L30;
  }
  csclr = *tol;
L30:
  c1r = cwr[1] * csclr;
  c1i = cwi[1] * csclr;
  c2r = cwr[2] * csclr;
  c2i = cwi[2] * csclr;
  str = yr[1];
  sti = yi[1];
  /* ----------------------------------------------------------------------- */
  /*     CINU=CINU*(CONJG(CT)/ABS(CT))*(1.0D0/ABS(CT) PREVENTS */
  /*     UNDER- OR OVERFLOW PREMATURELY BY SQUARING ABS(CT) */
  /* ----------------------------------------------------------------------- */
  ptr = str * c1r - sti * c1i;
  pti = str * c1i + sti * c1r;
  ptr += c2r;
  pti += c2i;
  ctr = *zrr * ptr - *zri * pti;
  cti = *zrr * pti + *zri * ptr;
  act = slatec_zabs_impl_(&ctr, &cti);
  ract = 1. / act;
  ctr *= ract;
  cti = -cti * ract;
  ptr = cinur * ract;
  pti = cinui * ract;
  cinur = ptr * ctr - pti * cti;
  cinui = ptr * cti + pti * ctr;
  yr[1] = cinur * csclr;
  yi[1] = cinui * csclr;
  if (*n == 1) {
    return 0;
  }
  i_u_1 = *n;
  for (i_u_ = 2; i_u_ <= i_u_1; ++i_u_) {
    ptr = str * cinur - sti * cinui;
    cinui = str * cinui + sti * cinur;
    cinur = ptr;
    str = yr[i_u_];
    sti = yi[i_u_];
    yr[i_u_] = cinur * csclr;
    yi[i_u_] = cinui * csclr;
    /* L40: */
  }
  return 0;
L50:
  *nz = -1;
  if (nw == -2) {
    *nz = -2;
  }
  return 0;
} /* slatec_zwrsk_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZWRSK_IMPL_H */