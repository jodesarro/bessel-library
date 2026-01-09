/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/slatec_zshch_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt

    Description:
        Slatec's ZSHCH subroutine, originally translated from Fortran to C
        with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZSHCH_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZSHCH_IMPL_H

#include <math.h>


/* zshch.f -- translated by f2c (version 20100827) */

/* DECK ZSHCH */
/* Subroutine */
static inline int slatec_zshch_impl_(double *zr, double *zi, double *cshr, 
	double *cshi, double *cchr, double *cchi)
{
    /* Local variables */
    static double ch, cn, sh, sn;

/* ***BEGIN PROLOGUE  ZSHCH */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESH and ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CSHCH-A, ZSHCH-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y) */
/*     AND CCH=COSH(X+I*Y), WHERE I**2=-1. */

/* ***SEE ALSO  ZBESH, ZBESK */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZSHCH */

/* ***FIRST EXECUTABLE STATEMENT  ZSHCH */
    sh = sinh(*zr);
    ch = cosh(*zr);
    sn = sin(*zi);
    cn = cos(*zi);
    *cshr = sh * cn;
    *cshi = ch * sn;
    *cchr = ch * cn;
    *cchi = sh * sn;
    return 0;
} /* slatec_zshch_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZSHCH_IMPL_H */