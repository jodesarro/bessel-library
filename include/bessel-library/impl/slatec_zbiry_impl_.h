/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/slatec_zbiry_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt

    Description:
        Slatec's ZBIRY subroutine, originally translated from Fortran to C
        with f2c, and adapted.
*/

#ifndef BESSEL_LIBRARY_SLATEC_ZBIRY_IMPL_H
#define BESSEL_LIBRARY_SLATEC_ZBIRY_IMPL_H

#include <math.h>
#include "slatec_zabs_impl_.h"
#include "fortran_d1mach_impl_.h"
#include "fortran_i1mach_impl_.h"
#include "slatec_zabs_impl_.h"
#include "slatec_zbinu_impl_.h"
#include "slatec_zdiv_impl_.h"
#include "slatec_zsqrt_impl_.h"
#include "fmax_impl_.h"
#include "fmin_impl_.h"
#include "min_impl_.h"
#include "f2c_pow_dd_impl_.h"

/* zbiry.f -- translated by f2c (version 20100827) */

/* DECK ZBIRY */
/* Subroutine */
static inline int slatec_zbiry_impl_(double *zr, double *zi, int *id, 
	int *kode, double *bir, double *bii, int *ierr)
{
    /* Initialized data */
    static int c_u_1 = 1;
    static int c_u_2 = 2;
    static int c_u_4 = 4;
    static int c_u_5 = 5;
    static int c_u_9 = 9;
    static int c_u_14 = 14;
    static int c_u_15 = 15;
    static int c_u_16 = 16;
    static double tth = .666666666666666667;
    static double c1 = .614926627446000736;
    static double c2 = .448288357353826359;
    static double coef = .577350269189625765;
    static double pi = 3.14159265358979324;
    static double coner = 1.;
    static double conei = 0.;

    /* System generated locals */
    int i_u_1, i_u_2;
    double d_u_1;

    /* Local variables */
    static int k;
    static double d1, d2;
    static int k1, k2;
    static double aa, bb, ad, cc, ak, bk, ck, dk, az, rl;
    static int nz;
    static double s1i, az3, s2i, s1r, s2r, z3i, z3r, eaa, fid, dig, cyi[3]
	    , fmr, r1m5, fnu, cyr[3], tol, sti, str, sfac, alim, elim, csqi;
    static double atrm, fnul, ztai, csqr;
    static double ztar, trm1i, trm2i, trm1r, trm2r;

/* ***BEGIN PROLOGUE  ZBIRY */
/* ***PURPOSE  Compute the Airy function Bi(z) or its derivative dBi/dz */
/*            for complex argument z.  A scaling option is available */
/*            to help avoid overflow. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C10D */
/* ***TYPE      COMPLEX (CBIRY-C, ZBIRY-C) */
/* ***KEYWORDS  AIRY FUNCTION, BESSEL FUNCTION OF ORDER ONE THIRD, */
/*             BESSEL FUNCTION OF ORDER TWO THIRDS */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*                      ***A DOUBLE PRECISION ROUTINE*** */
/*         On KODE=1, ZBIRY computes the complex Airy function Bi(z) */
/*         or its derivative dBi/dz on ID=0 or ID=1 respectively. */
/*         On KODE=2, a scaling option exp(abs(Re(zeta)))*Bi(z) or */
/*         exp(abs(Re(zeta)))*dBi/dz is provided to remove the */
/*         exponential behavior in both the left and right half planes */
/*         where zeta=(2/3)*z**(3/2). */

/*         The Airy functions Bi(z) and dBi/dz are analytic in the */
/*         whole z-plane, and the scaling option does not destroy this */
/*         property. */

/*         Input */
/*           ZR     - DOUBLE PRECISION real part of argument Z */
/*           ZI     - DOUBLE PRECISION imag part of argument Z */
/*           ID     - Order of derivative, ID=0 or ID=1 */
/*           KODE   - A parameter to indicate the scaling option */
/*                    KODE=1  returns */
/*                            BI=Bi(z)  on ID=0 */
/*                            BI=dBi/dz on ID=1 */
/*                            at z=Z */
/*                        =2  returns */
/*                            BI=exp(abs(Re(zeta)))*Bi(z)  on ID=0 */
/*                            BI=exp(abs(Re(zeta)))*dBi/dz on ID=1 */
/*                            at z=Z where zeta=(2/3)*z**(3/2) */

/*         Output */
/*           BIR    - DOUBLE PRECISION real part of result */
/*           BII    - DOUBLE PRECISION imag part of result */
/*           IERR   - Error flag */
/*                    IERR=0  Normal return     - COMPUTATION COMPLETED */
/*                    IERR=1  Input error       - NO COMPUTATION */
/*                    IERR=2  Overflow          - NO COMPUTATION */
/*                            (Re(Z) too large with KODE=1) */
/*                    IERR=3  Precision warning - COMPUTATION COMPLETED */
/*                            (Result has less than half precision) */
/*                    IERR=4  Precision error   - NO COMPUTATION */
/*                            (Result has no precision) */
/*                    IERR=5  Algorithmic error - NO COMPUTATION */
/*                            (Termination condition not met) */

/* *Long Description: */

/*         Bi(z) and dBi/dz are computed from I Bessel functions by */

/*                Bi(z) =  c*sqrt(z)*( I(-1/3,zeta) + I(1/3,zeta) ) */
/*               dBi/dz =  c*   z   *( I(-2/3,zeta) + I(2/3,zeta) ) */
/*                    c =  1/sqrt(3) */
/*                 zeta =  (2/3)*z**(3/2) */

/*         when abs(z)>1 and from power series when abs(z)<=1. */

/*         In most complex variable computation, one must evaluate ele- */
/*         mentary functions.  When the magnitude of Z is large, losses */
/*         of significance by argument reduction occur.  Consequently, if */
/*         the magnitude of ZETA=(2/3)*Z**(3/2) exceeds U1=SQRT(0.5/UR), */
/*         then losses exceeding half precision are likely and an error */
/*         flag IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is */
/*         double precision unit roundoff limited to 18 digits precision. */
/*         Also, if the magnitude of ZETA is larger than U2=0.5/UR, then */
/*         all significance is lost and IERR=4.  In order to use the INT */
/*         function, ZETA must be further restricted not to exceed */
/*         U3=I1MACH(9)=LARGEST INTEGER.  Thus, the magnitude of ZETA */
/*         must be restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, */
/*         and U3 are approximately 2.0E+3, 4.2E+6, 2.1E+9 in single */
/*         precision and 4.7E+7, 2.3E+15, 2.1E+9 in double precision. */
/*         This makes U2 limiting is single precision and U3 limiting */
/*         in double precision.  This means that the magnitude of Z */
/*         cannot exceed approximately 3.4E+4 in single precision and */
/*         2.1E+6 in double precision.  This also means that one can */
/*         expect to retain, in the worst cases on 32-bit machines, */
/*         no digits in single precision and only 6 digits in double */
/*         precision. */

/*         The approximate relative error in the magnitude of a complex */
/*         Bessel function can be expressed as P*10**S where P=MAX(UNIT */
/*         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre- */
/*         sents the increase in error due to argument reduction in the */
/*         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))), */
/*         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF */
/*         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may */
/*         have only absolute accuracy.  This is most likely to occur */
/*         when one component (in magnitude) is larger than the other by */
/*         several orders of magnitude.  If one component is 10**K larger */
/*         than the other, then one can expect only MAX(ABS(LOG10(P))-K, */
/*         0) significant digits; or, stated another way, when K exceeds */
/*         the exponent of P, no significant digits remain in the smaller */
/*         component.  However, the phase angle retains absolute accuracy */
/*         because, in complex arithmetic with precision P, the smaller */
/*         component will not (as a rule) decrease below P times the */
/*         magnitude of the larger component. In these extreme cases, */
/*         the principal phase angle is on the order of +P, -P, PI/2-P, */
/*         or -PI/2+P. */

/* ***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe- */
/*                 matical Functions, National Bureau of Standards */
/*                 Applied Mathematics Series 55, U. S. Department */
/*                 of Commerce, Tenth Printing (1972) or later. */
/*               2. D. E. Amos, Computation of Bessel Functions of */
/*                 Complex Argument and Large Order, Report SAND83-0643, */
/*                 Sandia National Laboratories, Albuquerque, NM, May */
/*                 1983. */
/*               3. D. E. Amos, A Subroutine Package for Bessel Functions */
/*                 of a Complex Argument and Nonnegative Order, Report */
/*                 SAND85-1018, Sandia National Laboratory, Albuquerque, */
/*                 NM, May 1985. */
/*               4. D. E. Amos, A portable package for Bessel functions */
/*                 of a complex argument and nonnegative order, ACM */
/*                 Transactions on Mathematical Software, 12 (September */
/*                 1986), pp. 265-273. */

/* ***ROUTINES CALLED  D1MACH, I1MACH, ZABS, ZBINU, ZDIV, ZSQRT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   890801  REVISION DATE from Version 3.2 */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/*   920128  Category corrected.  (WRB) */
/*   920811  Prologue revised.  (DWL) */
/*   930122  Added ZSQRT to EXTERNAL statement.  (RWC) */
/* ***END PROLOGUE  ZBIRY */
/*     COMPLEX BI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3 */
/* ***FIRST EXECUTABLE STATEMENT  ZBIRY */
    *ierr = 0;
    nz = 0;
    if (*id < 0 || *id > 1) {
	*ierr = 1;
    }
    if (*kode < 1 || *kode > 2) {
	*ierr = 1;
    }
    if (*ierr != 0) {
	return 0;
    }
    az = slatec_zabs_impl_(zr, zi);
/* Computing MAX */
    d_u_1 = fortran_d1mach_impl_(&c_u_4);
    tol = fmax_impl_(d_u_1,1e-18);
    fid = (double) (*id);
    if (az > 1.f) {
	goto L70;
    }
/* ----------------------------------------------------------------------- */
/*     POWER SERIES FOR ABS(Z).LE.1. */
/* ----------------------------------------------------------------------- */
    s1r = coner;
    s1i = conei;
    s2r = coner;
    s2i = conei;
    if (az < tol) {
	goto L130;
    }
    aa = az * az;
    if (aa < tol / az) {
	goto L40;
    }
    trm1r = coner;
    trm1i = conei;
    trm2r = coner;
    trm2i = conei;
    atrm = 1.;
    str = *zr * *zr - *zi * *zi;
    sti = *zr * *zi + *zi * *zr;
    z3r = str * *zr - sti * *zi;
    z3i = str * *zi + sti * *zr;
    az3 = az * aa;
    ak = fid + 2.;
    bk = 3. - fid - fid;
    ck = 4. - fid;
    dk = fid + 3. + fid;
    d1 = ak * dk;
    d2 = bk * ck;
    ad = fmin_impl_(d1,d2);
    ak = fid * 9. + 24.;
    bk = 30. - fid * 9.;
    for (k = 1; k <= 25; ++k) {
	str = (trm1r * z3r - trm1i * z3i) / d1;
	trm1i = (trm1r * z3i + trm1i * z3r) / d1;
	trm1r = str;
	s1r += trm1r;
	s1i += trm1i;
	str = (trm2r * z3r - trm2i * z3i) / d2;
	trm2i = (trm2r * z3i + trm2i * z3r) / d2;
	trm2r = str;
	s2r += trm2r;
	s2i += trm2i;
	atrm = atrm * az3 / ad;
	d1 += ak;
	d2 += bk;
	ad = fmin_impl_(d1,d2);
	if (atrm < tol * ad) {
	    goto L40;
	}
	ak += 18.;
	bk += 18.;
/* L30: */
    }
L40:
    if (*id == 1) {
	goto L50;
    }
    *bir = c1 * s1r + c2 * (*zr * s2r - *zi * s2i);
    *bii = c1 * s1i + c2 * (*zr * s2i + *zi * s2r);
    if (*kode == 1) {
	return 0;
    }
    slatec_zsqrt_impl_(zr, zi, &str, &sti);
    ztar = tth * (*zr * str - *zi * sti);
    ztai = tth * (*zr * sti + *zi * str);
    aa = ztar;
    aa = -fabs(aa);
    eaa = exp(aa);
    *bir *= eaa;
    *bii *= eaa;
    return 0;
L50:
    *bir = s2r * c2;
    *bii = s2i * c2;
    if (az <= tol) {
	goto L60;
    }
    cc = c1 / (fid + 1.);
    str = s1r * *zr - s1i * *zi;
    sti = s1r * *zi + s1i * *zr;
    *bir += cc * (str * *zr - sti * *zi);
    *bii += cc * (str * *zi + sti * *zr);
L60:
    if (*kode == 1) {
	return 0;
    }
    slatec_zsqrt_impl_(zr, zi, &str, &sti);
    ztar = tth * (*zr * str - *zi * sti);
    ztai = tth * (*zr * sti + *zi * str);
    aa = ztar;
    aa = -fabs(aa);
    eaa = exp(aa);
    *bir *= eaa;
    *bii *= eaa;
    return 0;
/* ----------------------------------------------------------------------- */
/*     CASE FOR ABS(Z).GT.1.0 */
/* ----------------------------------------------------------------------- */
L70:
    fnu = (fid + 1.) / 3.;
/* ----------------------------------------------------------------------- */
/*     SET PARAMETERS RELATED TO MACHINE CONSTANTS. */
/*     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18. */
/*     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT. */
/*     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND */
/*     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR */
/*     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE. */
/*     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z. */
/*     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG). */
/*     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU. */
/* ----------------------------------------------------------------------- */
    k1 = fortran_i1mach_impl_(&c_u_15);
    k2 = fortran_i1mach_impl_(&c_u_16);
    r1m5 = fortran_d1mach_impl_(&c_u_5);
/* Computing MIN */
    i_u_1 = abs(k1), i_u_2 = abs(k2);
    k = min_impl_(i_u_1,i_u_2);
    elim = (k * r1m5 - 3.) * 2.303;
    k1 = fortran_i1mach_impl_(&c_u_14) - 1;
    aa = r1m5 * k1;
    dig = fmin_impl_(aa,18.);
    aa *= 2.303;
/* Computing MAX */
    d_u_1 = -aa;
    alim = elim + fmax_impl_(d_u_1,-41.45);
    rl = dig * 1.2 + 3.;
    fnul = (dig - 3.) * 6. + 10.;
/* ----------------------------------------------------------------------- */
/*     TEST FOR RANGE */
/* ----------------------------------------------------------------------- */
    aa = .5 / tol;
    bb = fortran_i1mach_impl_(&c_u_9) * .5;
    aa = fmin_impl_(aa,bb);
    aa = f2c_pow_dd_impl_(&aa, &tth);
    if (az > aa) {
	goto L260;
    }
    aa = sqrt(aa);
    if (az > aa) {
	*ierr = 3;
    }
    slatec_zsqrt_impl_(zr, zi, &csqr, &csqi);
    ztar = tth * (*zr * csqr - *zi * csqi);
    ztai = tth * (*zr * csqi + *zi * csqr);
/* ----------------------------------------------------------------------- */
/*     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL */
/* ----------------------------------------------------------------------- */
    sfac = 1.;
    ak = ztai;
    if (*zr >= 0.) {
	goto L80;
    }
    bk = ztar;
    ck = -fabs(bk);
    ztar = ck;
    ztai = ak;
L80:
    if (*zi != 0. || *zr > 0.) {
	goto L90;
    }
    ztar = 0.;
    ztai = ak;
L90:
    aa = ztar;
    if (*kode == 2) {
	goto L100;
    }
/* ----------------------------------------------------------------------- */
/*     OVERFLOW TEST */
/* ----------------------------------------------------------------------- */
    bb = fabs(aa);
    if (bb < alim) {
	goto L100;
    }
    bb += log(az) * .25;
    sfac = tol;
    if (bb > elim) {
	goto L190;
    }
L100:
    fmr = 0.;
    if (aa >= 0. && *zr > 0.) {
	goto L110;
    }
    fmr = pi;
    if (*zi < 0.) {
	fmr = -pi;
    }
    ztar = -ztar;
    ztai = -ztai;
L110:
/* ----------------------------------------------------------------------- */
/*     AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA) */
/*     KODE=2 RETURNS EXP(-ABS(XZTA))*I(FNU,ZTA) FROM CBESI */
/* ----------------------------------------------------------------------- */
    slatec_zbinu_impl_(&ztar, &ztai, &fnu, kode, &c_u_1, &cyr[1], &cyi[1], &nz, &rl, &fnul, &tol, &
	    elim, &alim);
    if (nz < 0) {
	goto L200;
    }
    aa = fmr * fnu;
    z3r = sfac;
    str = cos(aa);
    sti = sin(aa);
    s1r = (str * cyr[1] - sti * cyi[1]) * z3r;
    s1i = (str * cyi[1] + sti * cyr[1]) * z3r;
    fnu = (2. - fid) / 3.;
    slatec_zbinu_impl_(&ztar, &ztai, &fnu, kode, &c_u_2, &cyr[1], &cyi[1], &nz, &rl, &fnul, &tol, &
	    elim, &alim);
    cyr[1] *= z3r;
    cyi[1] *= z3r;
    cyr[2] *= z3r;
    cyi[2] *= z3r;
/* ----------------------------------------------------------------------- */
/*     BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3 */
/* ----------------------------------------------------------------------- */
    slatec_zdiv_impl_(&cyr[1], &cyi[1], &ztar, &ztai, &str, &sti);
    s2r = (fnu + fnu) * str + cyr[2];
    s2i = (fnu + fnu) * sti + cyi[2];
    aa = fmr * (fnu - 1.);
    str = cos(aa);
    sti = sin(aa);
    s1r = coef * (s1r + s2r * str - s2i * sti);
    s1i = coef * (s1i + s2r * sti + s2i * str);
    if (*id == 1) {
	goto L120;
    }
    str = csqr * s1r - csqi * s1i;
    s1i = csqr * s1i + csqi * s1r;
    s1r = str;
    *bir = s1r / sfac;
    *bii = s1i / sfac;
    return 0;
L120:
    str = *zr * s1r - *zi * s1i;
    s1i = *zr * s1i + *zi * s1r;
    s1r = str;
    *bir = s1r / sfac;
    *bii = s1i / sfac;
    return 0;
L130:
    aa = c1 * (1. - fid) + fid * c2;
    *bir = aa;
    *bii = 0.;
    return 0;
L190:
    *ierr = 2;
    nz = 0;
    return 0;
L200:
    if (nz == -1) {
	goto L190;
    }
    nz = 0;
    *ierr = 5;
    return 0;
L260:
    *ierr = 4;
    nz = 0;
    return 0;
} /* slatec_zbiry_impl_ */

#endif /* BESSEL_LIBRARY_SLATEC_ZBIRY_IMPL_H */