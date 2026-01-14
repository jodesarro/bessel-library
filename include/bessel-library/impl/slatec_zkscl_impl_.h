/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/slatec_zkscl_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Slatec's ZKSCL subroutine, originally translated from Fortran to C
        with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZKSCL_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZKSCL_IMPL_H

#include <math.h>
#include "slatec_zabs_impl_.h"
#include "slatec_zlog_impl_.h"
#include "slatec_zuchk_impl_.h"
#include "min_impl_.h"

/* zkscl.f -- translated by f2c (version 20100827) */

/* DECK ZKSCL */
/* Subroutine */
static inline int slatec_zkscl_impl_(double *zrr, double *zri, double *fnu,
	 int *n, double *yr, double *yi, int *nz, double *
	rzr, double *rzi, double *ascle, double *tol, double *
	elim)
{
    /* Initialized data */
    static double zeror = 0.;
    static double zeroi = 0.;

    /* System generated locals */
    int i_u_1;

    /* Local variables */
    static int i_u_, ic;
    static double as, fn;
    static int kk, nn, nw;
    static double s1i, s2i, s1r, s2r, acs, cki, elm, csi, ckr, cyi[2], 
	    zdi, csr, cyr[2], zdr, str, alas;
    static int idum;
    static double helim, celmr;

/* ***BEGIN PROLOGUE  ZKSCL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CKSCL-A, ZKSCL-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE */
/*     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN */
/*     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL. */

/* ***SEE ALSO  ZBESK */
/* ***ROUTINES CALLED  ZABS, ZLOG, ZUCHK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/*   930122  Added ZLOG to EXTERNAL statement.  (RWC) */
/* ***END PROLOGUE  ZKSCL */
/*     COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  ZKSCL */
    *nz = 0;
    ic = 0;
    nn = min_impl_(2,*n);
    i_u_1 = nn;
    for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
	s1r = yr[i_u_];
	s1i = yi[i_u_];
	cyr[i_u_ - 1] = s1r;
	cyi[i_u_ - 1] = s1i;
	as = slatec_zabs_impl_(&s1r, &s1i);
	acs = -(*zrr) + log(as);
	++(*nz);
	yr[i_u_] = zeror;
	yi[i_u_] = zeroi;
	if (acs < -(*elim)) {
	    goto L10;
	}
	slatec_zlog_impl_(&s1r, &s1i, &csr, &csi, &idum);
	csr -= *zrr;
	csi -= *zri;
	str = exp(csr) / *tol;
	csr = str * cos(csi);
	csi = str * sin(csi);
	slatec_zuchk_impl_(&csr, &csi, &nw, ascle, tol);
	if (nw != 0) {
	    goto L10;
	}
	yr[i_u_] = csr;
	yi[i_u_] = csi;
	ic = i_u_;
	--(*nz);
L10:
	;
    }
    if (*n == 1) {
	return 0;
    }
    if (ic > 1) {
	goto L20;
    }
    yr[1] = zeror;
    yi[1] = zeroi;
    *nz = 2;
L20:
    if (*n == 2) {
	return 0;
    }
    if (*nz == 0) {
	return 0;
    }
    fn = *fnu + 1.;
    ckr = fn * *rzr;
    cki = fn * *rzi;
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    helim = *elim * .5;
    elm = exp(-(*elim));
    celmr = elm;
    zdr = *zrr;
    zdi = *zri;

/*     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF */
/*     S2 GETS LARGER THAN EXP(ELIM/2) */

    i_u_1 = *n;
    for (i_u_ = 3; i_u_ <= i_u_1; ++i_u_) {
	kk = i_u_;
	csr = s2r;
	csi = s2i;
	s2r = ckr * csr - cki * csi + s1r;
	s2i = cki * csr + ckr * csi + s1i;
	s1r = csr;
	s1i = csi;
	ckr += *rzr;
	cki += *rzi;
	as = slatec_zabs_impl_(&s2r, &s2i);
	alas = log(as);
	acs = -zdr + alas;
	++(*nz);
	yr[i_u_] = zeror;
	yi[i_u_] = zeroi;
	if (acs < -(*elim)) {
	    goto L25;
	}
	slatec_zlog_impl_(&s2r, &s2i, &csr, &csi, &idum);
	csr -= zdr;
	csi -= zdi;
	str = exp(csr) / *tol;
	csr = str * cos(csi);
	csi = str * sin(csi);
	slatec_zuchk_impl_(&csr, &csi, &nw, ascle, tol);
	if (nw != 0) {
	    goto L25;
	}
	yr[i_u_] = csr;
	yi[i_u_] = csi;
	--(*nz);
	if (ic == kk - 1) {
	    goto L40;
	}
	ic = kk;
	goto L30;
L25:
	if (alas < helim) {
	    goto L30;
	}
	zdr -= *elim;
	s1r *= celmr;
	s1i *= celmr;
	s2r *= celmr;
	s2i *= celmr;
L30:
	;
    }
    *nz = *n;
    if (ic == *n) {
	*nz = *n - 1;
    }
    goto L45;
L40:
    *nz = kk - 2;
L45:
    i_u_1 = *nz;
    for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
	yr[i_u_] = zeror;
	yi[i_u_] = zeroi;
/* L50: */
    }
    return 0;
} /* slatec_zkscl_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZKSCL_IMPL_H */