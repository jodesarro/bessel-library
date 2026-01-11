/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/core/cyl_k.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 with guards for C++98 compatibility
    References: include/bessel-library/references.txt

    Description:
        Computes, in double complex type for C, or in std::complex<double>
        type for C++, modified cylindrical Bessel functions of the second
        kind, real order, and complex argument.
*/

#ifndef BESSEL_LIBRARY_CYL_K_H
#define BESSEL_LIBRARY_CYL_K_H

#ifdef __cplusplus

/* Includes, typedefs and/or macros for C++98 compatibility */
#include <complex> /* for complex numbers */
typedef std::complex<double> tpdcomplex_impl_;

#else

#include <complex.h> /* for complex numbers */
typedef double complex tpdcomplex_impl_;

#endif /* __cplusplus */

#include "../impl/cyl_k_full_seq_impl_.h"

/*
    Returns the modified cylindrical Bessel function of the second kind, real
    order nu, and complex argument z, i.e., K_nu(z).

    Parameters:
    - nu, real order of K_nu(z).
    - z, complex argument of K_nu(z).
*/
static inline tpdcomplex_impl_ cyl_k(double nu, tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ ck[1];
    
    /* Compute cyl_k_full_seq_impl_ */
    cyl_k_full_seq_impl_(nu, 1, z, ck, 0);

    /* Return */
    return ck[0];
}

/*
    Returns the scaled version of the modified cylindrical Bessel function of
    the second kind, real order nu, and complex argument z, i.e.,
    K_nu(z)*exp(z).

    Parameters:
    - nu, real order of K_nu(z)*exp(z).
    - z, complex argument of K_nu(z)*exp(z).
*/
static inline tpdcomplex_impl_ cyl_k_scal(double nu,
    tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ ck[1];
    
    /* Compute cyl_k_full_seq_impl_ */
    cyl_k_full_seq_impl_(nu, 1, z, ck, 1);

    /* Return */
    return ck[0];
}

/*
    Computes a n-sequency array of modified cylindrical Bessel functions of
    the second kind, real order nu, and complex argument z, i.e.,
    {K_nu(z), Y_(nu+1)(z), ..., Y_(nu+n-1)(z)}.

    Parameters:
    - nu, real order of K_nu(z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_k_arr array.
    - z, complex argument of K_nu(z).
    - cyl_k_arr, array of size n to output K_nu(z) for the orders nu,
    nu+1, ..., nu+n-1
*/
static inline void cyl_k_seq(double nu, int n, tpdcomplex_impl_ z,
    tpdcomplex_impl_ *cyl_k_arr) {

    cyl_k_full_seq_impl_(nu, n, z, cyl_k_arr, 0);
}

/*
    Computes a n-sequency array of the scaled version of the modified
    cylindrical Bessel functions of the second kind, real order nu, and
    complex argument z, i.e., {K_nu(z)*exp(z), Y_(nu+1)(z)*exp(z),
    ..., Y_(nu+n-1)(z)*exp(z)}.

    Parameters:
    - nu, real order of K_nu(z)*exp(z).
    - n, number of elements in the sequence for computing the orders nu,
    nu+1, ..., nu+n-1 It is also the size of the cyl_k_scaled_arr array.
    - z, complex argument of K_nu(z)*exp(z).
    - cyl_k_scaled_arr, array of size n to output K_nu(z)*exp(z) for the
    orders nu, nu+1, ..., nu+n-1
*/
static inline void cyl_k_scal_seq(double nu, int n,
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_k_scaled_arr) {
    
    cyl_k_full_seq_impl_(nu, n, z, cyl_k_scaled_arr, 1);
}

#endif /* BESSEL_LIBRARY_CYL_K_H */