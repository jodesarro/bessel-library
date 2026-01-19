/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/cyl_i_full_seq_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Implements and computes a n-sequency array in double complex type for
        C, or in std::complex<double> for C++, of modified cylindrical Bessel
        functions of the first kind, real order nu, and complex argument z,
        i.e., {I_nu(z), I_(nu+1)(z), ..., I_(nu+n-1)(z)}, for also negative
        orders, by means of the routines from the Slatec library and
        recurrence relations for negative orders.
*/

#ifndef BESSEL_LIBRARY_CYL_I_FULL_SEQ_IMPL_H
#define BESSEL_LIBRARY_CYL_I_FULL_SEQ_IMPL_H

#include <stdlib.h> /* For malloc and free */
#include <math.h> /* For math operations */
#include <float.h>  /* For DBL_EPSILON */
#include "cplx_c_cpp_impl_.h"
#include "slatec_zbesi_impl_.h"
#include "slatec_zbesk_impl_.h"
#include "slatec_flags_impl_.h"

/* Fallback for M_PI if not defined by <math.h> */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
    Implements and computes a n-sequency array of modified cylindrical Bessel
    functions of the first kind, real order nu, and complex argument z, i.e.,
    {I_nu(z), I_(nu+1)(z), ..., I_(nu+n-1)(z)}, for also negative orders, by
    means of the routines from the Slatec library and recurrence relations for
    negative orders.
    
    Parameters:
    - nu, real order of I_nu(z).
    - n, number n of elements in the sequence for computing the orders nu,
    nu+1, ..., nu+n-1. It is also the size of the cyl_i_arr array.
    - z, complex argument of I_nu(z).
    - cyl_i_arr, array of size n to output I_nu(z) for the orders nu, nu+1,
    ..., nu+n-1
    - scaled, returns the scaled version I_nu(z)*exp(-abs(real(z))) if 1.

    Implementation: In general, the implementation is based on the D. E. Amos
    Fortran 77 routines from the Slatec library [3]. Such Fortran routines,
    and all their dependencies, were carefully translated to C. Negative
    orders are handled by Eqs. (6.1.5) and (6.5.4) of Ref. [2]
    for, respectively, nu integer and nu real; in
    the latter case, it yields INFINITY + I * INFINITY abs(z)=0.
*/
static inline void cyl_i_full_seq_impl_(double nu, int n,
    tpdfcplx_impl_ z, tpdfcplx_impl_ *cyl_i_arr, int scaled) {
    
    int kode = (scaled == 1 ? 2 : 1);
    int nz, ierr;
    double x = creal(z);
    double y = cimag(z);
    double fnu = fabs(nu);
    double nu_m = nu + (double)(n - 1);

    if (nu >= 0.0) {

        /* Positive orders */
        
        /* Dynamic mem alloc of auxiliary arrays */
        double *cir = (double*)malloc((n + 1) * sizeof(double));
        double *cii = (double*)malloc((n + 1) * sizeof(double));
        ++cir; ++cii;

        /* Compute zbesi */
        slatec_zbesi_impl_(&x, &y, &fnu, &kode, &n, &cir[0], &cii[0], &nz,
            &ierr);
        slatec_flags_zbesi_impl_(ierr, nz);
        
        /* Store in the array */
        for (int i = 0; i < n; i++) {
            cyl_i_arr[i] = cir[i] + I_IMPL_ * cii[i];
        }

        /* Free auxiliary pointers */
        free(--cir); free(--cii);
    }
    else if (nu_m <= 0.0) {
    
        /* Array of only negative orders */

        if (fabs(floor(fnu) - fnu) < DBL_EPSILON) {

            /* Integer negative orders */

            /* Dynamic mem alloc of auxiliary arrays */
            double *cir = (double*)malloc((n + 1) * sizeof(double));
            double *cii = (double*)malloc((n + 1) * sizeof(double));
            ++cir; ++cii;
            
            /* Compute zbesi */
            double fnu_m = fabs(nu_m);
            slatec_zbesi_impl_(&x, &y, &fnu_m, &kode, &n, &cir[0], &cii[0],
                &nz, &ierr);
            slatec_flags_zbesi_impl_(ierr, nz);
            
            /* Store in the array */
            for (int i = 0; i < n; i++) {
                int tmp = n - 1 - i;
                /* Eq. (6.1.5) of Ref. [2] */
                cyl_i_arr[i] = cir[tmp] + I_IMPL_ * cii[tmp];
            }

            /* Free auxiliary pointers */
            free(--cir); free(--cii);
        }
        else if (fabs(floor(nu_m) - nu_m) >= DBL_EPSILON
            && cabs(z) < DBL_EPSILON) {
        
            /* Non-int negative order with z=0 -> infinity */
        
            for (int i = 0; i < n; i++) {
                cyl_i_arr[i] = CPLX_IMPL_(INFINITY, INFINITY);
            }
        }
        else {
            
            /* Non-int negative orders with z!=0 */
            
            /* Dynamic mem alloc of auxiliary arrays */
            double *cir = (double*)malloc((n + 1) * sizeof(double));
            double *cii = (double*)malloc((n + 1) * sizeof(double));
            ++cir; ++cii;

            /* Compute zbesi */
            double fnu_m = fabs(nu_m);
            slatec_zbesi_impl_(&x, &y, &fnu_m, &kode, &n, &cir[0], &cii[0],
                &nz, &ierr);
            slatec_flags_zbesi_impl_(ierr, nz);

            /* Dynamic mem alloc of auxiliary arrays */
            double *ckr = (double*)malloc((n + 1) * sizeof(double));
            double *cki = (double*)malloc((n + 1) * sizeof(double));
            ++ckr; ++cki;

            /* Ensure non-scaled version for K_nu(z) */
            int kodek = 1;

            /* Compute zbesk */
            slatec_zbesk_impl_(&x, &y, &fnu_m, &kodek, &n, &ckr[0], &cki[0],
                &nz, &ierr);
            slatec_flags_zbesk_impl_(ierr, nz);

            /* Treat the scaled version of K_nu(z) for I_nu(z) */
            tpdfcplx_impl_ c_2_pi_e_m_abs_zr = 2.0/M_PI * (scaled == 1 ?
                                            exp(-fabs(creal(z))) : 1.0 );      

            /* Store in the array */
            for (int i = 0; i < n; i++) {
                int tmp1 = n - 1 - i;
                /* Eq. (6.5.4) of Ref. [2] */
                double tmp2 = (fnu_m + (double)tmp1) * M_PI;
                cyl_i_arr[i] = (cir[tmp1] + I_IMPL_ * cii[tmp1])
                         + (ckr[tmp1] + I_IMPL_ * cki[tmp1])
                         * sin(tmp2) * c_2_pi_e_m_abs_zr;
            }

            /* Free auxiliary pointers */
            free(--cir); free(--cii);
            free(--ckr); free(--cki);
        }
    }
    else {
        
        /* Negative and positive orders */

        int n_m = (int)floor(fabs(nu)) + 1;
        int n_p = n - n_m;
        double fnu_m = fabs(nu + (double)(n_m - 1));

        /* Negative orders */

        if (fabs(floor(fnu) - fnu) < DBL_EPSILON) {
            
            /* Integer negative orders */
            
            /* Dynamic mem alloc of auxiliary arrays */
            double *cir_m = (double*)malloc((n_m + 1) * sizeof(double));
            double *cii_m = (double*)malloc((n_m + 1) * sizeof(double));
            ++cir_m; ++cii_m;

            /* Compute zbesi */
            slatec_zbesi_impl_(&x, &y, &fnu_m, &kode, &n_m, &cir_m[0],
                &cii_m[0], &nz, &ierr);
            slatec_flags_zbesi_impl_(ierr, nz);

            /* Store in the array */
            for (int i = 0; i < n_m; i++) {
                int tmp = n_m - 1 - i;
                /* Eq. (6.1.5) of Ref. [2] */
                cyl_i_arr[i] = cir_m[tmp] + I_IMPL_ * cii_m[tmp];
            }

            /* Free auxiliary pointers */
            free(--cir_m); free(--cii_m);
        }
        else if (fabs(floor(nu_m) - nu_m) >= DBL_EPSILON
            && cabs(z) == 0.0) {
            
            /* Non-int negative order with z=0 -> infinity */
            
            for (int i = 0; i < n_m; i++) {
                cyl_i_arr[i] = CPLX_IMPL_(INFINITY, INFINITY);
            }
        }
        else {

            /* Non-int negative orders with z!=0 */

            /* Dynamic mem alloc of auxiliary arrays */
            double *cir_m = (double*)malloc((n_m + 1) * sizeof(double));
            double *cii_m = (double*)malloc((n_m + 1) * sizeof(double));
            ++cir_m; ++cii_m;
            
            /* Compute zbesi */
            slatec_zbesi_impl_(&x, &y, &fnu_m, &kode, &n_m, &cir_m[0],
                &cii_m[0], &nz, &ierr);
            slatec_flags_zbesi_impl_(ierr, nz);

            /* Dynamic mem alloc of auxiliary arrays */
            double *ckr_m = (double*)malloc((n_m + 1) * sizeof(double));
            double *cki_m = (double*)malloc((n_m + 1) * sizeof(double));
            ++ckr_m; ++cki_m;

            /* Ensure non-scaled version for K_nu(z) */
            int kodek = 1;

            /* Compute zbesk */
            slatec_zbesk_impl_(&x, &y, &fnu_m, &kodek, &n_m, &ckr_m[0],
                &cki_m[0], &nz, &ierr);
            slatec_flags_zbesk_impl_(ierr, nz);

            /* Treat the scaled version of K_nu(z) for I_nu(z) */
            tpdfcplx_impl_ c_2_pi_e_m_abs_zr = 2.0/M_PI * (scaled == 1 ?
                                            exp(-fabs(creal(z))) : 1.0 );     

            /* Store in the array */
            for (int i = 0; i < n_m; i++) {
                int tmp1 = n_m - 1 - i;
                /* Eq. (6.5.4) of Ref. [2] */
                double tmp2 = (fnu_m + (double)tmp1) * M_PI;
                cyl_i_arr[i] = (cir_m[tmp1] + I_IMPL_ * cii_m[tmp1])
                         + (ckr_m[tmp1] + I_IMPL_ * cki_m[tmp1])
                         * sin(tmp2) * c_2_pi_e_m_abs_zr;
            }

            /* Free auxiliary pointers */
            free(--cir_m); free(--cii_m);
            free(--ckr_m); free(--cki_m);
        }

        /* Positive orders */

        /* Dynamic mem alloc of auxiliary arrays */
        double *cir_p = (double*)malloc((n_p + 1) * sizeof(double));
        double *cii_p = (double*)malloc((n_p + 1) * sizeof(double));
        ++cir_p; ++cii_p;

        /* Compute zbesi */
        double fnu_p = nu + (double)n_m;
        slatec_zbesi_impl_(&x, &y, &fnu_p, &kode, &n_p, &cir_p[0], &cii_p[0],
            &nz, &ierr);
        slatec_flags_zbesi_impl_(ierr, nz);
        
        /* Store in the array */
        for (int i = n_m; i < n; i++) {
            int tmp1 = i - n_m;
            cyl_i_arr[i] = cir_p[tmp1] + I_IMPL_ * cii_p[tmp1];
        }

        /* Free auxiliary pointers */
        free(--cir_p); free(--cii_p);
    }
}

#endif /* BESSEL_LIBRARY_CYL_I_FULL_SEQ_IMPL_H */