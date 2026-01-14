/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/slatec_zdiv_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Slatec's ZDIV subroutine, originally translated from Fortran to C
        with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZDIV_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZDIV_IMPL_H

#include "slatec_zabs_impl_.h"

/* zdiv.f -- translated by f2c (version 20100827) */

/* DECK ZDIV */
/* Subroutine */
static inline int slatec_zdiv_impl_(double *ar, double *ai, double *br, 
	double *bi, double *cr, double *ci)
{
    static double ca, cb, cc, cd, bm;

/* ***BEGIN PROLOGUE  ZDIV */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and */
/*            ZBIRY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (ZDIV-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     DOUBLE PRECISION COMPLEX DIVIDE C=A/B. */

/* ***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY */
/* ***ROUTINES CALLED  ZABS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZDIV */
/* ***FIRST EXECUTABLE STATEMENT  ZDIV */
    bm = 1. / slatec_zabs_impl_(br, bi);
    cc = *br * bm;
    cd = *bi * bm;
    ca = (*ar * cc + *ai * cd) * bm;
    cb = (*ai * cc - *ar * cd) * bm;
    *cr = ca;
    *ci = cb;
    return 0;
} /* slatec_zdiv_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZDIV_IMPL_H */