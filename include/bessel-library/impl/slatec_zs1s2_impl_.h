/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/slatec_zs1s2_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt

    Description:
        Slatec's ZS1S2 subroutine, originally translated from Fortran to C
        with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZS1S2_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZS1S2_IMPL_H

#include <math.h>
#include "slatec_zabs_impl_.h"
#include "slatec_zexp_impl_.h"
#include "slatec_zlog_impl_.h"
#include "fmax_impl_.h"

/* zs1s2.f -- translated by f2c (version 20100827) */

/* DECK ZS1S2 */
/* Subroutine */
static inline int slatec_zs1s2_impl_(double *zrr, double *zri, double *s1r,
	 double *s1i, double *s2r, double *s2i, int *nz, 
	double *ascle, double *alim, int *iuf)
{
    /* Initialized data */
    static double zeror = 0.;
    static double zeroi = 0.;

    /* Local variables */
    static double aa, c1i, as1, as2, c1r, aln, s1di, s1dr;
    static int idum;

/* ***BEGIN PROLOGUE  ZS1S2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZAIRY and ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CS1S2-A, ZS1S2-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE */
/*     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON- */
/*     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION. */
/*     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF */
/*     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER */
/*     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE */
/*     PRECISION ABOVE THE UNDERFLOW LIMIT. */

/* ***SEE ALSO  ZAIRY, ZBESK */
/* ***ROUTINES CALLED  ZABS, ZEXP, ZLOG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/*   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC) */
/* ***END PROLOGUE  ZS1S2 */
/*     COMPLEX CZERO,C1,S1,S1D,S2,ZR */
/* ***FIRST EXECUTABLE STATEMENT  ZS1S2 */
    *nz = 0;
    as1 = slatec_zabs_impl_(s1r, s1i);
    as2 = slatec_zabs_impl_(s2r, s2i);
    if (*s1r == 0. && *s1i == 0.) {
	goto L10;
    }
    if (as1 == 0.) {
	goto L10;
    }
    aln = -(*zrr) - *zrr + log(as1);
    s1dr = *s1r;
    s1di = *s1i;
    *s1r = zeror;
    *s1i = zeroi;
    as1 = zeror;
    if (aln < -(*alim)) {
	goto L10;
    }
    slatec_zlog_impl_(&s1dr, &s1di, &c1r, &c1i, &idum);
    c1r = c1r - *zrr - *zrr;
    c1i = c1i - *zri - *zri;
    slatec_zexp_impl_(&c1r, &c1i, s1r, s1i);
    as1 = slatec_zabs_impl_(s1r, s1i);
    ++(*iuf);
L10:
    aa = fmax_impl_(as1,as2);
    if (aa > *ascle) {
	return 0;
    }
    *s1r = zeror;
    *s1i = zeroi;
    *s2r = zeror;
    *s2i = zeroi;
    *nz = 1;
    *iuf = 0;
    return 0;
} /* slatec_zs1s2_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZS1S2_IMPL_H */