/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/slatec_zuni2_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Slatec's ZUNI2 subroutine translated from Fortran to C with f2c.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZUNI2_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZUNI2_IMPL_H

#include <math.h>
#include "fortran_d1mach_impl_.h"
#include "slatec_zabs_impl_.h"
#include "slatec_zairy_impl_.h"
#include "slatec_zuchk_impl_.h"
#include "slatec_zunhj_impl_.h"
#include "slatec_zuoik_impl_.h"
#include "min_impl_.h"
#include "fmax_impl_.h"

/* zuni2.f -- translated by f2c (version 20100827) */

/* DECK ZUNI2 */
/* Subroutine */
static inline int slatec_zuni2_impl_(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *yr, double *yi, int *
	nz, int *nlast, double *fnul, double *tol, double *
	elim, double *alim)
{
    /* Initialized data */
    static int c_u_0 = 0;
    static int c_u_1 = 1;
    static int c_u_2 = 2;
    static double zeror = 0.;
    static double zeroi = 0.;
    static double coner = 1.;
    static double cipr[4] = { 1.,0.,-1.,0. };
    static double cipi[4] = { 0.,1.,0.,-1. };
    static double hpi = 1.57079632679489662;
    static double aic = 1.265512123484645396;

    /* System generated locals */
    int i_u_1;

    /* Local variables */
    static int i_u_, j, k, nd;
    static double fn;
    static int in, nn, nw;
    static double c2i, c2m, c1r, c2r, s1i, s2i, rs1, s1r, s2r, aii, ang, 
	    car;
    static int nai;
    static double air, zbi, cyi[2], sar;
    static int nuf, inu;
    static double bry[3], raz, sti, zbr, zni, cyr[2], rzi, str, znr, rzr, 
	    daii, cidi, aarg;
    static int ndai;
    static double dair, aphi, argi, cscl, phii, crsc, argr;
    static int idum;
    static double phir, csrr[3], cssr[3], rast;
    static int iflag;
    static double ascle, asumi, bsumi;
    static double asumr, bsumr;
    static double zeta1i, zeta2i, zeta1r, zeta2r;

/* ***BEGIN PROLOGUE  ZUNI2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to ZBESI and ZBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CUNI2-A, ZUNI2-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF */
/*     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I */
/*     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO. */

/*     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC */
/*     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET. */
/*     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER */
/*     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL. */
/*     Y(I)=CZERO FOR I=NLAST+1,N */

/* ***SEE ALSO  ZBESI, ZBESK */
/* ***ROUTINES CALLED  D1MACH, ZABS, ZAIRY, ZUCHK, ZUNHJ, ZUOIK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ZUNI2 */
/*     COMPLEX AI,ARG,ASUM,BSUM,CFN,CI,CID,CIP,CONE,CRSC,CSCL,CSR,CSS, */
/*    *CZERO,C1,C2,DAI,PHI,RZ,S1,S2,Y,Z,ZB,ZETA1,ZETA2,ZN */
    /* Parameter adjustments */
    --yi;
    --yr;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  ZUNI2 */
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
/*     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI */
/* ----------------------------------------------------------------------- */
    znr = *zi;
    zni = -(*zr);
    zbr = *zr;
    zbi = *zi;
    cidi = -coner;
    inu = (int) (*fnu);
    ang = hpi * (*fnu - inu);
    c2r = cos(ang);
    c2i = sin(ang);
    car = c2r;
    sar = c2i;
    in = inu + *n - 1;
    in = in % 4 + 1;
    str = c2r * cipr[in - 1] - c2i * cipi[in - 1];
    c2i = c2r * cipi[in - 1] + c2i * cipr[in - 1];
    c2r = str;
    if (*zi > 0.) {
	goto L10;
    }
    znr = -znr;
    zbi = -zbi;
    cidi = -cidi;
    c2i = -c2i;
L10:
/* ----------------------------------------------------------------------- */
/*     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER */
/* ----------------------------------------------------------------------- */
    fn = fmax_impl_(*fnu,1.);
    slatec_zunhj_impl_(&znr, &zni, &fn, &c_u_1, tol, &phir, &phii, &argr, &argi, &zeta1r, &
	    zeta1i, &zeta2r, &zeta2i, &asumr, &asumi, &bsumr, &bsumi);
    if (*kode == 1) {
	goto L20;
    }
    str = zbr + zeta2r;
    sti = zbi + zeta2i;
    rast = fn / slatec_zabs_impl_(&str, &sti);
    str = str * rast * rast;
    sti = -sti * rast * rast;
    s1r = -zeta1r + str;
    s1i = -zeta1i + sti;
    goto L30;
L20:
    s1r = -zeta1r + zeta2r;
    s1i = -zeta1i + zeta2i;
L30:
    rs1 = s1r;
    if (fabs(rs1) > *elim) {
	goto L150;
    }
L40:
    nn = min_impl_(2,nd);
    i_u_1 = nn;
    for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
	fn = *fnu + (nd - i_u_);
	slatec_zunhj_impl_(&znr, &zni, &fn, &c_u_0, tol, &phir, &phii, &argr, &argi, &
		zeta1r, &zeta1i, &zeta2r, &zeta2i, &asumr, &asumi, &bsumr, &
		bsumi);
	if (*kode == 1) {
	    goto L50;
	}
	str = zbr + zeta2r;
	sti = zbi + zeta2i;
	rast = fn / slatec_zabs_impl_(&str, &sti);
	str = str * rast * rast;
	sti = -sti * rast * rast;
	s1r = -zeta1r + str;
	s1i = -zeta1i + sti + fabs(*zi);
	goto L60;
L50:
	s1r = -zeta1r + zeta2r;
	s1i = -zeta1i + zeta2i;
L60:
/* ----------------------------------------------------------------------- */
/*     TEST FOR UNDERFLOW AND OVERFLOW */
/* ----------------------------------------------------------------------- */
	rs1 = s1r;
	if (fabs(rs1) > *elim) {
	    goto L120;
	}
	if (i_u_ == 1) {
	    iflag = 2;
	}
	if (fabs(rs1) < *alim) {
	    goto L70;
	}
/* ----------------------------------------------------------------------- */
/*     REFINE  TEST AND SCALE */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
	aphi = slatec_zabs_impl_(&phir, &phii);
	aarg = slatec_zabs_impl_(&argr, &argi);
	rs1 = rs1 + log(aphi) - log(aarg) * .25 - aic;
	if (fabs(rs1) > *elim) {
	    goto L120;
	}
	if (i_u_ == 1) {
	    iflag = 1;
	}
	if (rs1 < 0.) {
	    goto L70;
	}
	if (i_u_ == 1) {
	    iflag = 3;
	}
L70:
/* ----------------------------------------------------------------------- */
/*     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR */
/*     EXPONENT EXTREMES */
/* ----------------------------------------------------------------------- */
	slatec_zairy_impl_(&argr, &argi, &c_u_0, &c_u_2, &air, &aii, &nai, &idum);
	slatec_zairy_impl_(&argr, &argi, &c_u_1, &c_u_2, &dair, &daii, &ndai, &idum);
	str = dair * bsumr - daii * bsumi;
	sti = dair * bsumi + daii * bsumr;
	str += air * asumr - aii * asumi;
	sti += air * asumi + aii * asumr;
	s2r = phir * str - phii * sti;
	s2i = phir * sti + phii * str;
	str = exp(s1r) * cssr[iflag - 1];
	s1r = str * cos(s1i);
	s1i = str * sin(s1i);
	str = s2r * s1r - s2i * s1i;
	s2i = s2r * s1i + s2i * s1r;
	s2r = str;
	if (iflag != 1) {
	    goto L80;
	}
	slatec_zuchk_impl_(&s2r, &s2i, &nw, bry, tol);
	if (nw != 0) {
	    goto L120;
	}
L80:
	if (*zi <= 0.) {
	    s2i = -s2i;
	}
	str = s2r * c2r - s2i * c2i;
	s2i = s2r * c2i + s2i * c2r;
	s2r = str;
	cyr[i_u_ - 1] = s2r;
	cyi[i_u_ - 1] = s2i;
	j = nd - i_u_ + 1;
	yr[j] = s2r * csrr[iflag - 1];
	yi[j] = s2i * csrr[iflag - 1];
	str = -c2i * cidi;
	c2i = c2r * cidi;
	c2r = str;
/* L90: */
    }
    if (nd <= 2) {
	goto L110;
    }
    raz = 1. / slatec_zabs_impl_(zr, zi);
    str = *zr * raz;
    sti = -(*zi) * raz;
    rzr = (str + str) * raz;
    rzi = (sti + sti) * raz;
    bry[1] = 1. / bry[0];
    bry[2] = fortran_d1mach_impl_(&c_u_2);
    s1r = cyr[0];
    s1i = cyi[0];
    s2r = cyr[1];
    s2i = cyi[1];
    c1r = csrr[iflag - 1];
    ascle = bry[iflag - 1];
    k = nd - 2;
    fn = (double) k;
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
	    goto L100;
	}
	str = fabs(c2r);
	sti = fabs(c2i);
	c2m = fmax_impl_(str,sti);
	if (c2m <= ascle) {
	    goto L100;
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
L100:
	;
    }
L110:
    return 0;
L120:
    if (rs1 > 0.) {
	goto L140;
    }
/* ----------------------------------------------------------------------- */
/*     SET UNDERFLOW AND UPDATE PARAMETERS */
/* ----------------------------------------------------------------------- */
    yr[nd] = zeror;
    yi[nd] = zeroi;
    ++(*nz);
    --nd;
    if (nd == 0) {
	goto L110;
    }
    slatec_zuoik_impl_(zr, zi, fnu, kode, &c_u_1, &nd, &yr[1], &yi[1], &nuf, tol, elim, 
	    alim);
    if (nuf < 0) {
	goto L140;
    }
    nd -= nuf;
    *nz += nuf;
    if (nd == 0) {
	goto L110;
    }
    fn = *fnu + (nd - 1);
    if (fn < *fnul) {
	goto L130;
    }
/*      FN = CIDI */
/*      J = NUF + 1 */
/*      K = MOD(J,4) + 1 */
/*      S1R = CIPR(K) */
/*      S1I = CIPI(K) */
/*      IF (FN.LT.0.0D0) S1I = -S1I */
/*      STR = C2R*S1R - C2I*S1I */
/*      C2I = C2R*S1I + C2I*S1R */
/*      C2R = STR */
    in = inu + nd - 1;
    in = in % 4 + 1;
    c2r = car * cipr[in - 1] - sar * cipi[in - 1];
    c2i = car * cipi[in - 1] + sar * cipr[in - 1];
    if (*zi <= 0.) {
	c2i = -c2i;
    }
    goto L40;
L130:
    *nlast = nd;
    return 0;
L140:
    *nz = -1;
    return 0;
L150:
    if (rs1 > 0.) {
	goto L140;
    }
    *nz = *n;
    i_u_1 = *n;
    for (i_u_ = 1; i_u_ <= i_u_1; ++i_u_) {
	yr[i_u_] = zeror;
	yi[i_u_] = zeroi;
/* L160: */
    }
    return 0;
} /* slatec_zuni2_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZUNI2_IMPL_H */