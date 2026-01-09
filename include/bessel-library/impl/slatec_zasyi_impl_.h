/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/slatec_zasyi_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt

    Description:
        Slatec's ZASYI subroutine, originally translated from Fortran to C
        with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZASYI_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZASYI_IMPL_H

#include <math.h>
#include "fortran_d1mach_impl_.h"
#include "slatec_zabs_impl_.h"
#include "slatec_zdiv_impl_.h"
#include "slatec_zexp_impl_.h"
#include "slatec_zmlt_impl_.h"
#include "slatec_zsqrt_impl_.h"
#include "min_impl_.h"

/* DECK ZASYI */
/* Subroutine */
static inline int slatec_zasyi_impl_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *yr, double *yi, int *
	nz, double *rl, double *tol, double *elim, double *
	alim)
{
    /* Initialized data */
    static int c_u_1 = 1;
    static double pi = 3.14159265358979324;
    static double rtpi = .159154943091895336;
    static double zeror = 0.;
    static double zeroi = 0.;
    static double coner = 1.;
    static double conei = 0.;

    /* System generated locals */
    int i_u_1, i_u_2;
    double d_u_1, d_u_2;

    /* Local variables */
    static int i_u_, j, k, m;
    static double s, aa, bb;
    static int ib;
    static double ak, bk;
    static int il, jl;
    static double az;
    static int nn;
    static double p1i, s2i, p1r, s2r, cki, dki, fdn, arg, aez, arm, ckr, 
	    dkr, czi, ezi, sgn;
    static int inu;
    static double raz, czr, ezr, sqk, sti, rzi, tzi, str, rzr, tzr, ak1i, 
	    ak1r, cs1i, cs2i, cs1r, cs2r, dnu2, rtr1, dfnu, atol;
    static int koded;

/* ***BEGIN PROLOGUE  ZASYI */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CASYI-A, ZASYI-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY */
/*     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE ABS(Z) IN THE */
/*     REGION ABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN. */
/*     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1. */

/* ***SEE ALSO  ZBESI, ZBESK */
/* ***ROUTINES CALLED  D1MACH, ZABS, ZDIV, ZEXP, ZMLT, ZSQRT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/*   930122  Added ZEXP and ZSQRT to EXTERNAL statement.  (RWC) */
/* ***END PROLOGUE  ZASYI */
/*     COMPLEX AK1,CK,CONE,CS1,CS2,CZ,CZERO,DK,EZ,P1,RZ,S2,Y,Z */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  ZASYI */
    *nz = 0;
    az = slatec_zabs_impl_(zr, zi);
    arm = fortran_d1mach_impl_(&c_u_1) * 1e3;
    rtr1 = sqrt(arm);
    il = min_impl_(2,*n);
    dfnu = *fnu + (*n - il);
/* ----------------------------------------------------------------------- */
/*     OVERFLOW TEST */
/* ----------------------------------------------------------------------- */
    raz = 1. / az;
    str = *zr * raz;
    sti = -(*zi) * raz;
    ak1r = rtpi * str * raz;
    ak1i = rtpi * sti * raz;
    slatec_zsqrt_impl_(&ak1r, &ak1i, &ak1r, &ak1i);
    czr = *zr;
    czi = *zi;
    if (*kode != 2) {
	goto L10;
    }
    czr = zeror;
    czi = *zi;
L10:
    if (fabs(czr) > *elim) {
	goto L100;
    }
    dnu2 = dfnu + dfnu;
    koded = 1;
    if (fabs(czr) > *alim && *n > 2) {
	goto L20;
    }
    koded = 0;
    slatec_zexp_impl_(&czr, &czi, &str, &sti);
    slatec_zmlt_impl_(&ak1r, &ak1i, &str, &sti, &ak1r, &ak1i);
L20:
    fdn = 0.;
    if (dnu2 > rtr1) {
	fdn = dnu2 * dnu2;
    }
    ezr = *zr * 8.;
    ezi = *zi * 8.;
/* ----------------------------------------------------------------------- */
/*     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE */
/*     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE */
/*     EXPANSION FOR THE IMAGINARY PART. */
/* ----------------------------------------------------------------------- */
    aez = az * 8.;
    s = *tol / aez;
    jl = (int) (*rl + *rl + 2);
    p1r = zeror;
    p1i = zeroi;
    if (*zi == 0.) {
	goto L30;
    }
/* ----------------------------------------------------------------------- */
/*     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF */
/*     SIGNIFICANCE WHEN FNU OR N IS LARGE */
/* ----------------------------------------------------------------------- */
    inu = (int) (*fnu);
    arg = (*fnu - inu) * pi;
    inu = inu + *n - il;
    ak = -sin(arg);
    bk = cos(arg);
    if (*zi < 0.) {
	bk = -bk;
    }
    p1r = ak;
    p1i = bk;
    if (inu % 2 == 0) {
	goto L30;
    }
    p1r = -p1r;
    p1i = -p1i;
L30:
    i_u_1 = il;
    for (k = 1; k <= i_u_1; ++k) {
	sqk = fdn - 1.;
	atol = s * fabs(sqk);
	sgn = 1.;
	cs1r = coner;
	cs1i = conei;
	cs2r = coner;
	cs2i = conei;
	ckr = coner;
	cki = conei;
	ak = 0.;
	aa = 1.;
	bb = aez;
	dkr = ezr;
	dki = ezi;
	i_u_2 = jl;
	for (j = 1; j <= i_u_2; ++j) {
	    slatec_zdiv_impl_(&ckr, &cki, &dkr, &dki, &str, &sti);
	    ckr = str * sqk;
	    cki = sti * sqk;
	    cs2r += ckr;
	    cs2i += cki;
	    sgn = -sgn;
	    cs1r += ckr * sgn;
	    cs1i += cki * sgn;
	    dkr += ezr;
	    dki += ezi;
	    aa = aa * fabs(sqk) / bb;
	    bb += aez;
	    ak += 8.;
	    sqk -= ak;
	    if (aa <= atol) {
		goto L50;
	    }
/* L40: */
	}
	goto L110;
L50:
	s2r = cs1r;
	s2i = cs1i;
	if (*zr + *zr >= *elim) {
	    goto L60;
	}
	tzr = *zr + *zr;
	tzi = *zi + *zi;
	d_u_1 = -tzr;
	d_u_2 = -tzi;
	slatec_zexp_impl_(&d_u_1, &d_u_2, &str, &sti);
	slatec_zmlt_impl_(&str, &sti, &p1r, &p1i, &str, &sti);
	slatec_zmlt_impl_(&str, &sti, &cs2r, &cs2i, &str, &sti);
	s2r += str;
	s2i += sti;
L60:
	fdn = fdn + dfnu * 8. + 4.;
	p1r = -p1r;
	p1i = -p1i;
	m = *n - il + k;
	yr[m] = s2r * ak1r - s2i * ak1i;
	yi[m] = s2r * ak1i + s2i * ak1r;
/* L70: */
    }
    if (*n <= 2) {
	return 0;
    }
    nn = *n;
    k = nn - 2;
    ak = (double) k;
    str = *zr * raz;
    sti = -(*zi) * raz;
    rzr = (str + str) * raz;
    rzi = (sti + sti) * raz;
    ib = 3;
    i_u_1 = nn;
    for (i_u_ = ib; i_u_ <= i_u_1; ++i_u_) {
	yr[k] = (ak + *fnu) * (rzr * yr[k + 1] - rzi * yi[k + 1]) + yr[k + 2];
	yi[k] = (ak + *fnu) * (rzr * yi[k + 1] + rzi * yr[k + 1]) + yi[k + 2];
	ak += -1.;
	--k;
/* L80: */
    }
    if (koded == 0) {
	return 0;
    }
    slatec_zexp_impl_(&czr, &czi, &ckr, &cki);
    i_u_1 = nn;
    for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
	str = yr[i_u_] * ckr - yi[i_u_] * cki;
	yi[i_u_] = yr[i_u_] * cki + yi[i_u_] * ckr;
	yr[i_u_] = str;
/* L90: */
    }
    return 0;
L100:
    *nz = -1;
    return 0;
L110:
    *nz = -2;
    return 0;
} /* slatec_zasyi_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZASYI_IMPL_H */