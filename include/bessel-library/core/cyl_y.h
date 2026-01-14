/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/core/cyl_y.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 with guards for C++98 compatibility
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Computes, in double complex type for C, or in std::complex<double>
        type for C++, cylindrical Bessel functions of the second kind, real
        order, and complex argument.
*/

#ifndef BESSEL_LIBRARY_CYL_Y_H
#define BESSEL_LIBRARY_CYL_Y_H

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

#include "../impl/cyl_y_full_seq_impl_.h"

/*
    Returns the cylindrical Bessel function of the second kind, real order nu,
    and complex argument z, i.e., Y_nu(z).

    Parameters:
    - nu, real order of Y_nu(z).
    - z, complex argument of Y_nu(z).
    
    Implementation:
    - Similar to the cyl_y_seq() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdcomplex_impl_ cyl_y(double nu, tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ cy[1];
    
    /* Compute cyl_y_full_seq_impl_ */
    cyl_y_full_seq_impl_(nu, 1, z, cy, 0);

    /* Return */
    return cy[0];
}

/*
    Returns the scaled version of the cylindrical Bessel function of the
    second kind, real order nu, and complex argument z, i.e.,
    Y_nu(z)*exp(-abs(imag(z))).

    Parameters:
    - nu, real order of Y_nu(z)*exp(-abs(imag(z))).
    - z, complex argument of Y_nu(z)*exp(-abs(imag(z))).
        
    Implementation:
    - Similar to the cyl_y_seq() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdcomplex_impl_ cyl_y_scal(double nu, tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ cy[1];
    
    /* Compute cyl_y_full_seq_impl_ */
    cyl_y_full_seq_impl_(nu, 1, z, cy, 1);

    /* Return */
    return cy[0];
}

/*
    Computes a n-sequency array of cylindrical Bessel functions of the second
    kind, real order nu, and complex argument z, i.e.,
    {Y_nu(z), Y_(nu+1)(z), ..., Y_(nu+n-1)(z)}.

    Parameters:
    - nu, real order of Y_nu(z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_y_arr array.
    - z, complex argument of Y_nu(z).
    - cyl_y_arr, array of size n to output Y_nu(z) for the orders nu,
    nu+1, ..., nu+n-1.
    
    Implementation:
    - In general, the implementation is based on the D. E. Amos Fortran 77
    routines of the Slatec library [3]. Such Fortran routines,
    and all their dependencies, were carefully translated to C. Negative
    orders are handled by Eqs. (5.4.2) and (5.5.4) of Ref. [2]
    for, respectively, nu integer and nu real.
    When abs(z)=0, it yields -INFINITY if nu=0, or INFINITY + I * INFINITY
    otherwise.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
void cyl_y_seq(double nu, int n, tpdcomplex_impl_ z,
    tpdcomplex_impl_ *cyl_y_arr) {

    cyl_y_full_seq_impl_(nu, n, z, cyl_y_arr, 0);
}

/*
    Computes a n-sequency array of the scaled version of the cylindrical
    Bessel functions of the second kind, real order nu, and complex argument
    z, i.e., {Y_nu(z)*exp(-abs(imag(z))), Y_(nu+1)(z)*exp(-abs(imag(z))),
    ..., Y_(nu+n-1)(z)*exp(-abs(imag(z)))}.

    Parameters:
    - nu, real order of Y_nu(z)*exp(-abs(imag(z))).
    - n, number of elements in the sequence for computing the orders nu,
    nu+1, ..., nu+n-1 It is also the size of the cyl_y_scaled_arr array.
    - z, complex argument of Y_nu(z)*exp(-abs(imag(z))).
    - cyl_y_scaled_arr, array of size n to output Y_nu(z)*exp(-abs(imag(z)))
    for the orders nu, nu+1, ..., nu+n-1.
        
    Implementation:
    - Similar to the cyl_y_seq() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
void cyl_y_scal_seq(double nu, int n,
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_y_scaled_arr) {
    
    cyl_y_full_seq_impl_(nu, n, z, cyl_y_scaled_arr, 1);
}

#endif /* BESSEL_LIBRARY_CYL_Y_H */