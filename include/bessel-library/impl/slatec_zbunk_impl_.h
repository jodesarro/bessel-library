/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/slatec_zbunk_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Slatec's ZBUNK subroutine, originally translated from Fortran to C
        with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZBUNK_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZBUNK_IMPL_H

#include <math.h>
#include "slatec_zunk1_impl_.h"
#include "slatec_zunk2_impl_.h"

/* zbunk.f -- translated by f2c (version 20100827) */

/* DECK ZBUNK */
/* Subroutine */
static inline int slatec_zbunk_impl_(double *zr, double *zi, double *fnu, 
	int *kode, int *mr, int *n, double *yr, double *
	yi, int *nz, double *tol, double *elim, double *alim)
{
    static double ax, ay;

/* ***BEGIN PROLOGUE  ZBUNK */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESH and ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CBUNI-A, ZBUNI-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL. */
/*     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z) */
/*     IN ZUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN ZUNK2 */

/* ***SEE ALSO  ZBESH, ZBESK */
/* ***ROUTINES CALLED  ZUNK1, ZUNK2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZBUNK */
/*     COMPLEX Y,Z */
/* ***FIRST EXECUTABLE STATEMENT  ZBUNK */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */
    *nz = 0;
    ax = fabs(*zr) * 1.7321;
    ay = fabs(*zi);
    if (ay > ax) {
	goto L10;
    }
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN */
/*     -PI/3.LE.ARG(Z).LE.PI/3 */
/* ----------------------------------------------------------------------- */
    slatec_zunk1_impl_(zr, zi, fnu, kode, mr, n, &yr[1], &yi[1], nz, tol, elim, alim);
    goto L20;
L10:
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU */
/*     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I */
/*     AND HPI=PI/2 */
/* ----------------------------------------------------------------------- */
    slatec_zunk2_impl_(zr, zi, fnu, kode, mr, n, &yr[1], &yi[1], nz, tol, elim, alim);
L20:
    return 0;
} /* slatec_zbunk_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZBUNK_IMPL_H */