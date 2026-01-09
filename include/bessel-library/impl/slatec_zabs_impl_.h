/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/slatec_zabs_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt

    Description:
        Slatec's ZABS subroutine, originally translated from Fortran to C
        with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZABS_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZABS_IMPL_H

#include <math.h>

/* zabs.f -- translated by f2c (version 20100827) */

/* DECK ZABS */
/* Subroutine */
static inline double slatec_zabs_impl_(double *zr, double *zi)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    static double q, s, u, v;

/* ***BEGIN PROLOGUE  ZABS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and */
/*            ZBIRY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (ZABS-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE */
/*     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI) */

/* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZABS */
/* ***FIRST EXECUTABLE STATEMENT  ZABS */
    u = fabs(*zr);
    v = fabs(*zi);
    s = u + v;
/* ----------------------------------------------------------------------- */
/*     S*1.0D0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A */
/*     TRUE FLOATING ZERO */
/* ----------------------------------------------------------------------- */
    s *= 1.;
    if (s == 0.) {
	goto L20;
    }
    if (u > v) {
	goto L10;
    }
    q = u / v;
    ret_val = v * sqrt(q * q + 1.);
    return ret_val;
L10:
    q = v / u;
    ret_val = u * sqrt(q * q + 1.);
    return ret_val;
L20:
    ret_val = 0.;
    return ret_val;
} /* slatec_zabs_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZABS_IMPL_H */