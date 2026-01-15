/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/core/cyl_h1.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 with guards for C++98 compatibility
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Computes, in double complex type for C, or in std::complex<double>
        type for C++, cylindrical Hankel functions of the first kind, real
        order, and complex argument.
*/

#ifndef BESSEL_LIBRARY_CYL_H1_H
#define BESSEL_LIBRARY_CYL_H1_H

#ifndef BESSEL_LIBRARY_STATIC_INLINE_IMPL_
#define BESSEL_LIBRARY_STATIC_INLINE_IMPL_ static inline
#endif

#ifdef __cplusplus

/* Includes, typedefs and/or macros for C++98 compatibility */

#include <complex> /* For complex numbers */
typedef std::complex<double> tpdcomplex_impl_;

#else

#include <complex.h> /* For complex numbers */
typedef double complex tpdcomplex_impl_;

#endif /* __cplusplus */

#include "../impl/cyl_h1_full_seq_impl_.h"

/*
    Returns the cylindrical Hankel function of the first kind, real order nu,
    and complex argument z, i.e., H1_nu(z).

    Parameters:
    - nu, real order of H1_nu(z).
    - z, complex argument of H1_nu(z).
            
    Implementation: Similar to the cyl_h1_seq() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdcomplex_impl_ cyl_h1(double nu, tpdcomplex_impl_ z) {

    /* Array of one size */
    tpdcomplex_impl_ ch1[1];
    
    /* Compute cyl_h1_full_seq_impl_ */
    cyl_h1_full_seq_impl_(nu, 1, z, ch1, 0);

    /* Return */
    return ch1[0];
}

/*
    Returns the scaled version of the cylindrical Hankel function of the
    first kind, real order nu, and complex argument z, i.e.,
    H1_nu(z)*exp(-i*z).

    Parameters:
    - nu, real order of H1_nu(z)*exp(-i*z).
    - z, complex argument of H1_nu(z)*exp(-i*z).
                
    Implementation: Similar to the cyl_h1_seq() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdcomplex_impl_ cyl_h1_scal(double nu, tpdcomplex_impl_ z) {

    /* Array of one size */
    tpdcomplex_impl_ ch1[1];
    
    /* Compute cyl_h1_full_seq_impl_ */
    cyl_h1_full_seq_impl_(nu, 1, z, ch1, 1);

    /* Return */
    return ch1[0];
}

/*
    Computes a n-sequency array of cylindrical Hankel functions of the first
    kind, real order nu, and complex argument z, i.e., {H1_nu(z),
    H1_(nu+1)(z), ..., H1_(nu+n-1)(z)}.

    Parameters:
    - nu, real order of H1_nu(z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_h1_arr array.
    - z, complex argument of H1_nu(z).
    - cyl_h1_arr, array of size n to output H1_nu(z) for the orders nu, nu+1,
    ..., nu+n-1.

    Implementation: In general, the implementation is based on the D. E. Amos
    Fortran 77 routines of the Slatec library [3]. Such Fortran routines,
    and all their dependencies, were carefully translated to C. Negative
    orders are handled by Eq. (9.1.6) of Ref. [1].
    It yields INFINITY + I * INFINITY when abs(z)=0.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
void cyl_h1_seq(double nu, int n,
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_h1_arr) {

    cyl_h1_full_seq_impl_(nu, n, z, cyl_h1_arr, 0);
}

/*
    Computes a n-sequency array of the scaled version of the cylindrical
    Hankel functions of the first kind, real order nu, and complex argument z,
    i.e., {H1_nu(z)*exp(-i*z), H1_(nu+1)(z)*exp(-i*z), ...,
    H1_(nu+n-1)(z)*exp(-i*z)}.

    Parameters:
    - nu, real order of H1_nu(z)*exp(-i*z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_h2_scal_arr array.
    - z, complex argument of H1_nu(z)*exp(-i*z).
    - cyl_h2_scal_arr, array of size n to output  H1_nu(z)*exp(-i*z) for the
    orders nu, nu+1, ..., nu+n-1.
                
    Implementation: Similar to the cyl_h1_seq() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
void cyl_h1_scal_seq(double nu, int n,
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_h2_scal_arr) {
    
    cyl_h1_full_seq_impl_(nu, n, z, cyl_h2_scal_arr, 1);
}

#endif /* BESSEL_LIBRARY_CYL_H1_H */