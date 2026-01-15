/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/cyl_h2_full_seq_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 with guards for C++98 compatibility
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Implements and computes a n-sequency array, in double complex type for
        C, or in std::complex<double> type for C++, of cylindrical Hankel
        functions of the second kind, real order nu, and complex argument z,
        i.e., {H2_nu(z), H2_(nu+1)(z), ..., H2_(nu+n-1)(z)}, for also negative
        orders, by means of the routines of the Slatec library and recurrence
        relations for negative orders.
*/

#ifndef BESSEL_LIBRARY_CYL_H2_FULL_SEQ_IMPL_H
#define BESSEL_LIBRARY_CYL_H2_FULL_SEQ_IMPL_H

#ifdef __cplusplus

/* Includes, typedefs and/or macros for C++98 compatibility */

#include <complex> /* For complex numbers */
#include <cstdlib> /* For malloc and free */
typedef std::complex<double> tpdcomplex_impl_;
#define CPLX_IMPL_(x, y) tpdcomplex_impl_(x, y)
#define I_IMPL_ std::complex<double>(0.0, 1.0)
#define creal(z) std::real(z)
#define cimag(z) std::imag(z)
#define cabs(z) std::abs(z)
#define cexp(z) std::exp(z)

extern "C" {

#else

#include <complex.h> /* For complex numbers */
#include <stdlib.h> /* For malloc and free */
typedef double complex tpdcomplex_impl_;
#define I_IMPL_ I
#define CPLX_IMPL_(x, y) (x + I * y)

#endif /* __cplusplus */

#include <math.h> /* for math operations */
#include <float.h>  /* for DBL_EPSILON */
#include "slatec_zbesh_impl_.h"
#include "slatec_flags_impl_.h"

/* Fallback for M_PI if not defined by <math.h> */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
    Implements and computes a n-sequency array of cylindrical Hankel functions
    of the second kind, real order nu, and complex argument z, i.e.,
    {H2_nu(z), H2_(nu+1)(z), ..., H2_(nu+n-1)(z)}, for also negative orders, by
    means of the routines of the Slatec library and recurrence relations for
    negative orders.
    
    Parameters:
    - nu, real order of H2_nu(z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_h2_arr array.
    - z, complex argument of H2_nu(z).
    - cyl_h2_arr, array of size n to output H2_nu(z) for the orders nu, nu+1,
    ..., nu+n-1
    - scaled, returns the scaled version H2_nu(z)*exp(i*z)) if 1.
                
    Implementation:
    - In general, the implementation is based on the D. E. Amos Fortran 77
    routines of the Slatec library [3]. Such Fortran routines,
    and all their dependencies, were carefully translated to C. Negative
    orders are handled by Eq. (9.1.6) of Ref. [1]. It yields
    INFINITY + I * INFINITY when abs(z)=0.
*/
static inline void cyl_h2_full_seq_impl_(double nu, int n, tpdcomplex_impl_ z,
    tpdcomplex_impl_ *cyl_h2_arr, int scaled) {

    int kode = (scaled == 1 ? 2 : 1);
    int nz, ierr;
    double x = creal(z);
    double y = cimag(z);
    double fnu = fabs(nu);
    double nu_m = nu + (double)(n - 1);
    int m = 2;

    if (cabs(z) < DBL_EPSILON) {

        /* Complex infinity for all orders */
        for (int i = 0; i < n; i++) {
            cyl_h2_arr[i] = CPLX_IMPL_(INFINITY, INFINITY);
        }
    }
    else if (nu >= 0.0) {

        /* Dynamic mem alloc of auxiliary arrays */
        double *ch2r = (double*)malloc((n + 1) * sizeof(double));
        double *ch2i = (double*)malloc((n + 1) * sizeof(double));
        ++ch2r; ++ch2i;

        /* Compute zbesh */
        slatec_zbesh_impl_(&x, &y, &fnu, &kode, &m, &n, &ch2r[0], &ch2i[0],
            &nz, &ierr);
        slatec_flags_zbesh_impl_(ierr, nz);
        
        /* Store in the array */
        for (int i = 0; i < n; i++) {
            cyl_h2_arr[i] = ch2r[i] + I_IMPL_ * ch2i[i];
        }

        /* Free auxiliary pointers */
        free(--ch2r); free(--ch2i);
    }
    else if (nu_m <= 0.0) {
        
        /* Array of only negative orders */

        /* Dynamic mem alloc of auxiliary arrays */
        double *ch2r_m = (double*)malloc((n + 1) * sizeof(double));
        double *ch2i_m = (double*)malloc((n + 1) * sizeof(double));
        ++ch2r_m; ++ch2i_m;

        /* Compute zbesh */
        double fnu_m = fabs(nu_m);
        slatec_zbesh_impl_(&x, &y, &fnu_m, &kode, &m, &n, &ch2r_m[0],
            &ch2i_m[0], &nz, &ierr);
        slatec_flags_zbesh_impl_(ierr, nz);
        
        /* Store in the array */
        for (int i = 0; i < n; i++) {
            int tmp1 = n - 1 - i;
            double tmp2 = (fnu_m + (double)tmp1) * M_PI;
            /* Eq. (9.1.6) of Ref. [1] */
            cyl_h2_arr[i] = (ch2r_m[tmp1] + I_IMPL_ * ch2i_m[tmp1])
                          * cexp(-I_IMPL_ * tmp2);
        }
    
        /* Free auxiliary pointers */
        free(--ch2r_m); free(--ch2i_m);
    }
    else {

        /* Array of negative and positive orders */
        
        int n_m = (int)floor(fabs(nu)) + 1;
        int n_p = n - n_m;
        double fnu_m = fabs(nu + (double)(n_m - 1));

        /* Negative orders */

        /* Dynamic mem alloc of auxiliary arrays */
        double *ch2r_m = (double*)malloc((n_m + 1) * sizeof(double));
        double *ch2i_m = (double*)malloc((n_m + 1) * sizeof(double));
        ++ch2r_m; ++ch2i_m;

        /* Compute zbesh */
        slatec_zbesh_impl_(&x, &y, &fnu_m, &kode, &m, &n_m, &ch2r_m[0],
            &ch2i_m[0], &nz, &ierr);
        slatec_flags_zbesh_impl_(ierr, nz);
        
        /* Store in the array */
        for (int i = 0; i < n_m; i++) {
            int tmp1 = n_m - 1 - i;
            double tmp2 = (fnu_m + (double)tmp1) * M_PI;
            /* Eq. (9.1.6) of Ref. [1] */
            cyl_h2_arr[i] = (ch2r_m[tmp1] + I_IMPL_ * ch2i_m[tmp1])
                          * cexp(-I_IMPL_ * tmp2);
        }

        /* Free auxiliary pointers */
        free(--ch2r_m); free(--ch2i_m);

        /* Positive orders */

        /* Dynamic mem alloc of auxiliary arrays */
        double *ch2r_p = (double*)malloc((n_p + 1) * sizeof(double));
        double *ch2i_p = (double*)malloc((n_p + 1) * sizeof(double));
        ++ch2r_p; ++ch2i_p;

        /* Compute zbesh */
        double fnu_p = nu + (double)n_m;
        slatec_zbesh_impl_(&x, &y, &fnu_p, &kode, &m, &n_p, &ch2r_p[0],
            &ch2i_p[0], &nz, &ierr);
        slatec_flags_zbesh_impl_(ierr, nz);

        /* Store in the array */
        for (int i = n_m; i < n; i++) {
            int tmp1 = i - n_m;
            cyl_h2_arr[i] = ch2r_p[tmp1] + I_IMPL_ * ch2i_p[tmp1];
        }

        /* Free auxiliary pointers */
        free(--ch2r_p); free(--ch2i_p);
    }
}

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* BESSEL_LIBRARY_CYL_H2_FULL_SEQ_IMPL_H */