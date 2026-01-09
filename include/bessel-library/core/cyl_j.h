/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/core/cyl_j.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 with guards for C++98 compatibility
    References: include/bessel-library/references.txt

    Description:
        Computes, in double complex type for C, or in std::complex<double>
        type for C++, cylindrical Bessel functions of the first kind, real
        order, and complex argument.
*/

#ifndef BESSEL_LIBRARY_CYL_J_H
#define BESSEL_LIBRARY_CYL_J_H

#ifdef __cplusplus

/* Includes, typedefs and/or macros for C++98 compatibility */
#include <complex> /* for complex numbers */
typedef std::complex<double> tpdcomplex_impl_;

#else

#include <complex.h> /* for complex numbers */
typedef double complex tpdcomplex_impl_;

#endif /* __cplusplus */

#include "../impl/cyl_j_full_seq_impl_.h"

/*
    Returns the cylindrical Bessel function of the first kind, real order nu,
    and complex argument z, i.e., J_nu(z).

    Parameters:
    - nu, real order of J_nu(z).
    - z, complex argument of J_nu(z).
*/
static inline tpdcomplex_impl_ cyl_j(double nu, tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ cj[1];
    
    /* Compute cyl_j_full_seq_impl_ */
    cyl_j_full_seq_impl_(nu, 1, z, cj, 0);

    /* Return */
    return cj[0];
}

/*
    Returns the scaled version of the cylindrical Bessel function of the
    first kind, real order nu, and complex argument z, i.e.,
    J_nu(z)*exp(-abs(imag(z))).

    Parameters:
    - nu, real order of J_nu(z)*exp(-abs(imag(z))).
    - z, complex argument of J_nu(z)*exp(-abs(imag(z))).
*/
static inline tpdcomplex_impl_ cyl_j_scal(double nu,
    tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ cj[1];
    
    /* Compute cyl_j_full_seq_impl_ */
    cyl_j_full_seq_impl_(nu, 1, z, cj, 1);

    /* Return */
    return cj[0];
}

/*
    Computes a n-sequency array of cylindrical Bessel functions of the first
    kind, real order nu, and complex argument z, i.e.,
    {J_nu(z), J_(nu+1)(z), ..., J_(nu+n-1)(z)}.

    Parameters:
    - nu, real order of J_nu(z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_j_arr array.
    - z, complex argument of J_nu(z).
    - cyl_j_arr, array of size n to output J_nu(z) for the orders nu, nu+1,
    ..., nu+n-1
*/
static inline void cyl_j_seq(double nu, unsigned int n, tpdcomplex_impl_ z,
    tpdcomplex_impl_ *cyl_j_arr) {

    cyl_j_full_seq_impl_(nu, n, z, cyl_j_arr, 0);
}

/*
    Computes a n-sequency array of the scaled version of the cylindrical
    Bessel functions of the first kind, real order nu, and complex argument z,
    i.e., {J_nu(z)*exp(-abs(imag(z))), J_(nu+1)(z)*exp(-abs(imag(z))), ...,
    J_(nu+n-1)(z)*exp(-abs(imag(z)))}.

    Parameters:
    - nu, real order of J_nu(z)*exp(-abs(imag(z))).
    - n, number of elements in the sequence for computing the orders nu,
    nu+1, ..., nu+n-1 It is also the size of the cyl_j_scal_arr array.
    - z, complex argument of J_nu(z)*exp(-abs(imag(z))).
    - cyl_j_scal_arr, array of size n to output J_nu(z)*exp(-abs(imag(z)))
    for the orders nu, nu+1, ..., nu+n-1
*/
static inline void cyl_j_scal_seq(double nu, unsigned int n,
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_j_scal_arr) {
    
    cyl_j_full_seq_impl_(nu, n, z, cyl_j_scal_arr, 1);
}

#endif /* BESSEL_LIBRARY_CYL_J_H */