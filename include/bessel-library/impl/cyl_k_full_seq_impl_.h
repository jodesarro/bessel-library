/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/cyl_k_full_seq_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 with guards for C++98 compatibility
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Implements and computes a n-sequency array in double complex type for
        C, or in std::complex<double> for C++, of modified cylindrical Bessel
        functions of the second kind, real order nu, and complex argument z,
        i.e., {K_nu(z), K_(nu+1)(z), ..., K_(nu+n-1)(z)}, for also negative
        orders, by means of the routines of the Slatec library and recurrence
        relations for negative orders.
*/

#ifndef BESSEL_LIBRARY_CYL_K_FULL_SEQ_IMPL_H
#define BESSEL_LIBRARY_CYL_K_FULL_SEQ_IMPL_H

#ifdef __cplusplus

/* Includes, typedefs and/or macros for C++98 compatibility */

#include <complex> /* For complex numbers */
#include <cstdlib> /* For malloc and free */
typedef std::complex<double> tpdcomplex_impl_;
#define CPLX_impl_(x, y) tpdcomplex_impl_(x, y)
#define I_impl_ std::complex<double>(0.0, 1.0)
#define creal(z) std::real(z)
#define cimag(z) std::imag(z)
#define cabs(z) std::abs(z)

extern "C" {

#else

#include <complex.h> /* For complex numbers */
#include <stdlib.h> /* For malloc and free */
typedef double complex tpdcomplex_impl_;
#define I_impl_ I
#define CPLX_impl_(x, y) (x + I * y)

#endif /* __cplusplus */

#include <math.h> /* for math operations */
#include <float.h>  /* for DBL_EPSILON */
#include "slatec_zbesk_impl_.h"
#include "slatec_flags_impl_.h"

/* Fallback for M_PI if not defined by <math.h> */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
    Implements and computes a n-sequency array of modified cylindrical Bessel
    functions of the second kind, real order nu, and complex argument z, i.e.,
    {K_nu(z), K_(nu+1)(z), ..., K_(nu+n-1)(z)}, for also negative orders, by
    means of the routines of the Slatec library and recurrence relations for
    negative orders.
    
    Parameters:
    - nu, real order of K_nu(z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_k_arr array.
    - z, complex argument of K_nu(z).
    - cyl_k_arr, array of size n to output K_nu(z) for the orders nu, nu+1,
    ..., nu+n-1
    - scaled, returns the scaled version K_nu(z)*exp(z) if 1.
    
    Implementation:
    - In general, the implementation is based on the D. E. Amos Fortran 77
    routines of the Slatec library [3]. Such Fortran routines,
    and all their dependencies, were carefully translated to C. Negative
    orders are handled by Eqs. (6.5.5) of Ref. [2]. When
    abs(z)=0, it yields INFINITY if nu=0, or INFINITY + I * INFINITY
    otherwise.
*/
static inline void cyl_k_full_seq_impl_(double nu, int n, tpdcomplex_impl_ z,
    tpdcomplex_impl_ *cyl_k_arr, int scaled) {
    
    int kode = (scaled == 1 ? 2 : 1);
    int nz, ierr;
    double x = creal(z);
    double y = cimag(z);
    double fnu = fabs(nu);
    double nu_m = nu + (double)(n - 1);

    if (cabs(z) < DBL_EPSILON) {
        
        /* Complex infinity for all orders */
        for (int i = 0; i < n; i++) {
            cyl_k_arr[i] = CPLX_impl_(INFINITY, INFINITY);
        }
        if (fabs(nu_m - floor(nu_m)) < DBL_EPSILON && n > fnu) {
            
            /* Infinity for integer order 0 */
            cyl_k_arr[(int)floor(fnu)] = CPLX_impl_(INFINITY, 0.0);
        }
    }
    else if (nu >= 0.0) {

        /* Positive orders */
        
        /* Dynamic mem alloc of auxiliary arrays */
        double *ckr = (double*)malloc((n + 1) * sizeof(double));
        double *cki = (double*)malloc((n + 1) * sizeof(double));
        ++ckr; ++cki;
        
        /* Compute zbesk */
        slatec_zbesk_impl_(&x, &y, &fnu, &kode, &n, &ckr[0], &cki[0], &nz,
            &ierr);
        slatec_flags_zbesk_impl_(ierr, nz);

        /* Store in the array */
        for (int i = 0; i < n; i++) {
            cyl_k_arr[i] = ckr[i] + I_impl_ * cki[i];
        }
        
        /* Free auxiliary pointers */
        free(--ckr); free(--cki);
    }
    else if (nu_m <= 0.0) {
        
        /* Only negative orders */
        
        double *ckr = (double*)malloc((n + 1) * sizeof(double));
        double *cki = (double*)malloc((n + 1) * sizeof(double));
        ++ckr; ++cki;
        
        /* Compute zbesk */
        double fnu_m = fabs(nu_m);
        slatec_zbesk_impl_(&x, &y, &fnu_m, &kode, &n, &ckr[0], &cki[0], &nz,
            &ierr);
        slatec_flags_zbesk_impl_(ierr, nz);

        /* Store in the array */
        for (int i = 0; i < n; i++) {
            int tmp = n - 1 - i;
            /* Eq. (6.5.5) of Ref. [2] */
            cyl_k_arr[i] = ckr[tmp] + I_impl_ * cki[tmp];
        }
        
        /* Free auxiliary pointers */        
        free(--ckr); free(--cki);
    }
    else {
    
        /* Negative and positive orders */
        
        int n_m = (int)floor(fabs(nu)) + 1;
        int n_p = n - n_m;
        double fnu_m = fabs(nu + (double)(n_m - 1));

        /* Negative orders */

        /* Dynamic mem alloc of auxiliary arrays */
        double *ckr_m = (double*)malloc((n_m + 1) * sizeof(double));
        double *cki_m = (double*)malloc((n_m + 1) * sizeof(double));
        ++ckr_m; ++cki_m;

        /* Compute zbesk */
        slatec_zbesk_impl_(&x, &y, &fnu_m, &kode, &n_m, &ckr_m[0], &cki_m[0],
            &nz, &ierr);
        slatec_flags_zbesk_impl_(ierr, nz);

        /* Store in the array */
        for (int i = 0; i < n_m; i++) {
            int tmp = n_m - 1 - i;
            /* Eq. (6.5.5) of Ref. [2] */
            cyl_k_arr[i] = ckr_m[tmp] + I_impl_ * cki_m[tmp];
        }

        /* Free auxiliary pointers */
        free(--ckr_m); free(--cki_m);

        /* Positive orders */

        /* Dynamic mem alloc of auxiliary arrays */
        double *ckr_p = (double*)malloc((n_p + 1) * sizeof(double));
        double *cki_p = (double*)malloc((n_p + 1) * sizeof(double));
        ++ckr_p; ++cki_p;

        /* Compute zbesk */
        double fnu_p = nu + (double)n_m;
        slatec_zbesk_impl_(&x, &y, &fnu_p, &kode, &n_p, &ckr_p[0], &cki_p[0],
            &nz, &ierr);
        slatec_flags_zbesk_impl_(ierr, nz);
        
        /* Store in the array */
        for (int i = n_m; i < n; i++) {
            int tmp1 = i - n_m;
            cyl_k_arr[i] = ckr_p[tmp1] + I_impl_ * cki_p[tmp1];
        }
        
        /* Free auxiliary pointers */
        free(--ckr_p); free(--cki_p);
    }
}

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* BESSEL_LIBRARY_CYL_K_FULL_SEQ_IMPL_H */