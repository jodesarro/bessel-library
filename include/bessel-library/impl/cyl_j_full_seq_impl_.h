/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/cyl_j_full_seq_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 with guards for C++98 compatibility
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Implements and computes a n-sequency array in double complex type for
        C, or in std::complex<double> for C++, of cylindrical Bessel functions
        of the first kind, real order nu, and complex argument z, i.e.,
        {J_nu(z), J_(nu+1)(z), ..., J_(nu+n-1)(z)}, for also negative orders,
        by means of the routines of the Slatec library and recurrence
        relations for negative orders.
*/

#ifndef BESSEL_LIBRARY_CYL_J_FULL_SEQ_IMPL_H
#define BESSEL_LIBRARY_CYL_J_FULL_SEQ_IMPL_H

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
#define sin(z) std::sin(z)
#define cos(z) std::cos(z)

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
#include "slatec_zbesj_impl_.h"
#include "slatec_zbesy_impl_.h"
#include "slatec_flags_impl_.h"

/* Fallback for M_PI if not defined by <math.h> */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
    Implements and computes a n-sequency array of cylindrical Bessel functions
    of the first kind, real order nu, and complex argument z, i.e.,
    {J_nu(z), J_(nu+1)(z), ..., J_(nu+n-1)(z)}, for also negative orders, by
    means of the routines of the Slatec library and recurrence relations for
    negative orders.
    
    Parameters:
    - nu, real order of J_nu(z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_j_arr array.
    - z, complex argument of J_nu(z).
    - cyl_j_arr, array of size n to output J_nu(z) for the orders nu, nu+1,
    ..., nu+n-1
    - scaled, returns the scaled version J_nu(z)*exp(-abs(imag(z))) if 1.

    Implementation:
    - In general, the implementation is based on the D. E. Amos Fortran 77
    routines of the Slatec library [3]. Such Fortran routines,
    and all their dependencies, were carefully translated to C. Negative
    orders are handled by Eqs. (5.4.2) and (5.5.4) of Ref. [2]
    for, respectively, nu integer and nu real; in
    the latter case, it yields INFINITY + I * INFINITY when abs(z)=0.
*/
static inline void cyl_j_full_seq_impl_(double nu, int n, 
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_j_arr, int scaled) {

    int kode = (scaled == 1 ? 2 : 1);
    int nz, ierr;
    double x = creal(z);
    double y = cimag(z);
    double fnu = fabs(nu);
    double nu_m = nu + (double)(n - 1);

    if (nu >= 0.0) {
        
        /* Positive orders */

        /* Dynamic mem alloc of auxiliary arrays */
        double *cjr_ptr = (double*)malloc((n + 1) * sizeof(double));
        double *cji_ptr = (double*)malloc((n + 1) * sizeof(double));
        ++cjr_ptr; ++cji_ptr;

        /* Compute zbesj */
        slatec_zbesj_impl_(&x, &y, &fnu, &kode, &n, &cjr_ptr[0], &cji_ptr[0],
            &nz, &ierr);
        slatec_flags_zbesj_impl_(ierr, nz);
        
        /* Store in the array */
        for (int i = 0; i < n; i++) {
            cyl_j_arr[i] = cjr_ptr[i] + I_impl_ * cji_ptr[i];
        }

        /* Free auxiliary pointers */
        free(--cjr_ptr); free(--cji_ptr);
    }
    else if (nu_m <= 0.0) {
        
        /* Only negative orders */
        if (fabs(floor(fnu) - fnu) < DBL_EPSILON) {
            
            /* Integer negative orders */
            
            /* Dynamic mem alloc of auxiliary arrays */
            double *cjr_ptr = (double*)malloc((n + 1) * sizeof(double));
            double *cji_ptr = (double*)malloc((n + 1) * sizeof(double));
            ++cjr_ptr; ++cji_ptr;

            /* Compute zbesj */
            double fnu_m = fabs(nu_m);
            slatec_zbesj_impl_(&x, &y, &fnu_m, &kode, &n, &cjr_ptr[0],
                &cji_ptr[0], &nz, &ierr);
            slatec_flags_zbesj_impl_(ierr, nz);
            
            /* Store in the array */
            for (int i = 0; i < n; i++) {
                int tmp = n - 1 - i;
                cyl_j_arr[i] = pow(-1.0, fnu_m + tmp)
                         * (cjr_ptr[tmp] + I_impl_ * cji_ptr[tmp]);
            }

            /* Free auxiliary pointers */
            free(--cjr_ptr); free(--cji_ptr);

        }
        else if (fabs(floor(nu_m) - nu_m) >= DBL_EPSILON
                 && cabs(z) < DBL_EPSILON) {
            
            /* Non-int negative orders with abs(z) = 0 */
            for (int i = 0; i < n; i++) {
                cyl_j_arr[i] = CPLX_impl_(INFINITY, INFINITY);
            }
        }
        else {

            /* Non-int negative orders with abs(z) != 0 */

            /* Dynamic mem alloc of auxiliary arrays */
            double *cjr_ptr = (double*)malloc((n + 1) * sizeof(double));
            double *cji_ptr = (double*)malloc((n + 1) * sizeof(double));
            ++cjr_ptr; ++cji_ptr;

            /* Compute zbesj */
            double fnu_m = fabs(nu_m);
            slatec_zbesj_impl_(&x, &y, &fnu_m, &kode, &n, &cjr_ptr[0],
                &cji_ptr[0], &nz, &ierr);
            slatec_flags_zbesj_impl_(ierr, nz);

            /* Dynamic mem alloc of auxiliary arrays */
            double *cyr_ptr = (double*)malloc((n + 1) * sizeof(double));
            double *cyi_ptr = (double*)malloc((n + 1) * sizeof(double));
            double *cwrkr_ptr = (double*)malloc((n + 1) * sizeof(double));
            double *cwrki_ptr = (double*)malloc((n + 1) * sizeof(double));
            ++cyr_ptr; ++cyi_ptr; ++cwrkr_ptr; ++cwrki_ptr;

            /* Compute zbesy */
            slatec_zbesy_impl_(&x, &y, &fnu_m, &kode, &n, &cyr_ptr[0],
                &cyi_ptr[0], &nz, &cwrkr_ptr[0], &cwrki_ptr[0], &ierr);
            slatec_flags_zbesy_impl_(ierr, nz);

            /* Store in the array */
            for (int i = 0; i < n; i++) {
                int tmp1 = n - 1 - i;
                double tmp2 = (fnu_m + (double)tmp1) * M_PI;
                /* Eq. (5.5.4) of Ref. [2] */
                cyl_j_arr[i] = (cjr_ptr[tmp1] + I_impl_ * cji_ptr[tmp1])
                         * cos(tmp2)
                         - (cyr_ptr[tmp1] + I_impl_ * cyi_ptr[tmp1])
                         * sin(tmp2);
            }

            /* Free auxiliary pointers */
            free(--cjr_ptr); free(--cji_ptr);
            free(--cyr_ptr); free(--cyi_ptr);
            free(--cwrkr_ptr); free(--cwrki_ptr);
        }
    }
    else {

        /* Mixed negative and positive orders */
        int n_m = (int)floor(fabs(nu)) + 1;
        int n_p = n - n_m;
        double fnu_m = fabs(nu + (double)(n_m - 1));

        /* Negative orders */
        if (fabs(floor(fnu) - fnu) < DBL_EPSILON) {
            
            /* Dynamic mem alloc of auxiliary arrays */
            double *cjr_m_ptr = (double*)malloc((n_m + 1) * sizeof(double));
            double *cji_m_ptr = (double*)malloc((n_m + 1) * sizeof(double));
            ++cjr_m_ptr; ++cji_m_ptr;
            
            /* Compute zbesj */
            slatec_zbesj_impl_(&x, &y, &fnu_m, &kode, &n_m, &cjr_m_ptr[0],
                &cji_m_ptr[0], &nz, &ierr);
            slatec_flags_zbesj_impl_(ierr, nz);
            
            /* Store in the array */
            for (int i = 0; i < n_m; i++) {
                int tmp = n_m - 1 - i;
                cyl_j_arr[i] = pow(-1.0, fnu_m + tmp)
                         * (cjr_m_ptr[tmp] + I_impl_ * cji_m_ptr[tmp]);
            }
            
            /* Free auxiliary pointers */
            free(--cjr_m_ptr); free(--cji_m_ptr);

        }
        else if (fabs(floor(nu_m) - nu_m) >= DBL_EPSILON
                 && cabs(z) < DBL_EPSILON) {
            
            /* Store in the array */
            for (int i = 0; i < n_m; i++) {
                cyl_j_arr[i] = CPLX_impl_(INFINITY, INFINITY);
            }
        }
        else {

            /* Dynamic mem alloc of auxiliary arrays */
            double *cjr_m_ptr = (double*)malloc((n_m + 1) * sizeof(double));
            double *cji_m_ptr = (double*)malloc((n_m + 1) * sizeof(double));
            ++cjr_m_ptr; ++cji_m_ptr;

            /* Compute zbesj */
            slatec_zbesj_impl_(&x, &y, &fnu_m, &kode, &n_m, &cjr_m_ptr[0],
                &cji_m_ptr[0], &nz, &ierr);
            slatec_flags_zbesj_impl_(ierr, nz);

            /* Dynamic mem alloc of auxiliary arrays */
            double *cyr_m_ptr = (double*)malloc((n_m + 1) * sizeof(double));
            double *cyi_m_ptr = (double*)malloc((n_m + 1) * sizeof(double));
            double *cwrkr_ptr = (double*)malloc((n_m + 1) * sizeof(double));
            double *cwrki_ptr = (double*)malloc((n_m + 1) * sizeof(double));
            ++cyr_m_ptr; ++cyi_m_ptr; ++cwrkr_ptr; ++cwrki_ptr;

            /* Compute zbesy */
            slatec_zbesy_impl_(&x, &y, &fnu_m, &kode, &n_m, &cyr_m_ptr[0],
                &cyi_m_ptr[0], &nz, &cwrkr_ptr[0], &cwrki_ptr[0], &ierr);
            slatec_flags_zbesy_impl_(ierr, nz);

            /* Store in the array */
            for (int i = 0; i < n_m; i++) {
                int tmp1 = n_m - 1 - i;
                double tmp2 = (fnu_m + (double)tmp1) * M_PI;
                /* Eq. (5.5.4) of Ref. [2] */
                cyl_j_arr[i] = (cjr_m_ptr[tmp1] + I_impl_ * cji_m_ptr[tmp1])
                         * cos(tmp2)
                         - (cyr_m_ptr[tmp1] + I_impl_ * cyi_m_ptr[tmp1])
                         * sin(tmp2);
            }

            /* Free auxiliary pointers */
            free(--cjr_m_ptr); free(--cji_m_ptr);
            free(--cyr_m_ptr); free(--cyi_m_ptr);
            free(--cwrkr_ptr); free(--cwrki_ptr);
        }

        /* Positive orders */

        /* Dynamic mem alloc of auxiliary arrays */
        double *cjr_p_ptr = (double*)malloc((n_p + 1) * sizeof(double));
        double *cji_p_ptr = (double*)malloc((n_p + 1) * sizeof(double));
        ++cjr_p_ptr; ++cji_p_ptr;

        /* Compute zbesj */
        double fnu_p = nu + (double)n_m;
        slatec_zbesj_impl_(&x, &y, &fnu_p, &kode, &n_p, &cjr_p_ptr[0],
            &cji_p_ptr[0], &nz, &ierr);
        slatec_flags_zbesj_impl_(ierr, nz);

        /* Store in the array */
        for (int i = n_m; i < n; i++) {
            int tmp1 = i - n_m;
            cyl_j_arr[i] = cjr_p_ptr[tmp1] + I_impl_ * cji_p_ptr[tmp1];
        }

        /* Free auxiliary pointers */
        free(--cjr_p_ptr); free(--cji_p_ptr);
    }
}

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* BESSEL_LIBRARY_CYL_J_FULL_SEQ_IMPL_H */