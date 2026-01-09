/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/airy_ai_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 with guards for C++98 compatibility
    References: include/bessel-library/references.txt

    Description:
        Returns the Airy function of the first kind and complex argument z,
        i.e., Ai(z), in double complex type for C, or in std::complex<double> 
        type for C++, by means of the routines of the Slatec library.
*/

#ifndef BESSEL_LIBRARY_AIRY_AI_IMPL_H
#define BESSEL_LIBRARY_AIRY_AI_IMPL_H

#ifdef __cplusplus

/* Includes, typedefs and/or macros for C++98 compatibility */
#include <complex> /* for complex numbers */
typedef std::complex<double> tpdcomplex_impl_;
#define I_impl_ std::complex<double>(0.0, 1.0)
#define creal(z) std::real(z)
#define cimag(z) std::imag(z)

extern "C" {

#else

#include <complex.h> /* for complex numbers */
typedef double complex tpdcomplex_impl_;
#define I_impl_ I

#endif /* __cplusplus */

#include "slatec_zairy_impl_.h"
#include "slatec_flags_impl_.h"

/*
    Returns the Airy function of the first kind and complex argument z,
    i.e., Ai(z) by means of the routines of the Slatec library.
    
    Parameters:
    - z, complex argument of Ai(z).
    - derivative, replaces Ai(z) by dAi(z)/dz if 1.
    - scaled, returns the scaled version Ai(z)*exp((2/3)*pow(z,3/2)) if 1.
*/
static inline tpdcomplex_impl_ airy_ai_impl_(tpdcomplex_impl_ z,
    int derivative, int scaled) {
    
    int kode = ( scaled == 1 ? 2 : 1 );
    int nz, ierr;
    double x = creal(z);
    double y = cimag(z);

    /* Auxiliary arrays */
    double air_arr[2] = {0.0, 0.0}, aii_arr[2] = {0.0, 0.0};
    
    /* Compute zairy */
    slatec_zairy_impl_(&x, &y, &derivative, &kode, &air_arr[1], &aii_arr[1],
        &nz, &ierr);
    slatec_flags_zairy_impl_(ierr, nz);

    /* Return */
    return air_arr[1] + I_impl_ * aii_arr[1];
}

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* BESSEL_LIBRARY_AIRY_AI_IMPL_H */