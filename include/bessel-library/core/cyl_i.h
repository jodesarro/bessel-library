/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/core/cyl_i.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 with guards for C++98 compatibility
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Computes, in double complex type for C, or in std::complex<double>
        type for C++, modified cylindrical Bessel functions of the first kind,
        real order, and complex argument.
*/

#ifndef BESSEL_LIBRARY_CYL_I_H
#define BESSEL_LIBRARY_CYL_I_H

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

#include "../impl/cyl_i_full_seq_impl_.h"

/*
    Returns the modified cylindrical Bessel function of the first kind, real
    order nu, and complex argument z, i.e., I_nu(z).

    Parameters:
    - nu, real order of I_nu(z).
    - z, complex argument of I_nu(z).

    Implementation: Similar to the cyl_i_seq() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdcomplex_impl_ cyl_i(double nu, tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ ci[1];
    
    /* Compute cyl_i_full_seq_impl_ */
    cyl_i_full_seq_impl_(nu, 1, z, ci, 0);

    /* Return */
    return ci[0];
}

/*
    Returns the scaled version of the modified cylindrical Bessel function of
    the first kind, real order nu, and complex argument z, i.e.,
    I_nu(z)*exp(-abs(real(z))).

    Parameters:
    - nu, real order of I_nu(z)*exp(-abs(real(z))).
    - z, complex argument of I_nu(z)*exp(-abs(real(z))).
    
    Implementation: Similar to the cyl_i_seq() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdcomplex_impl_ cyl_i_scal(double nu, tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ ci[1];
    
    /* Compute cyl_i_full_seq_impl_ */
    cyl_i_full_seq_impl_(nu, 1, z, ci, 1);

    /* Return */
    return ci[0];
}

/*
    Computes a n-sequency array of modified cylindrical Bessel functions of
    the first kind, real order nu, and complex argument z, i.e.,
    {I_nu(z), I_(nu+1)(z), ..., I_(nu+n-1)(z)}.

    Parameters:
    - nu, real order of I_nu(z).
    - n, number n of elements in the sequence for computing the orders nu,
    nu+1, ..., nu+n-1. It is also the size of the cyl_i_arr array.
    - z, complex argument of I_nu(z).
    - cyl_i_arr, array of size n to output I_nu(z) for the orders nu, nu+1,
    ..., nu+n-1.

    Implementation: In general, the implementation is based on the D. E. Amos
    Fortran 77 routines from the Slatec library [3]. Such Fortran routines,
    and all their dependencies, were carefully translated to C. Negative
    orders are handled by Eqs. (6.1.5) and (6.5.4) of Ref. [2]
    for, respectively, nu integer and nu real; in
    the latter case, it yields INFINITY + I * INFINITY abs(z)=0.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
void cyl_i_seq(double nu, int n, tpdcomplex_impl_ z,
    tpdcomplex_impl_ *cyl_i_arr) {

    cyl_i_full_seq_impl_(nu, n, z, cyl_i_arr, 0);
}

/*
    Computes a n-sequency array of scaled versions of modified
    cylindrical Bessel functions of the first kind, real order nu, and complex
    argument z, i.e., {I_nu(z)*exp(-abs(real(z))),
    I_(nu+1)(z)*exp(-abs(real(z))), ..., I_(nu+n-1)(z)*exp(-abs(real(z)))}.

    Parameters:
    - nu, real order of I_nu(z)*exp(-abs(real(z))).
    - n, number n of elements in the sequence for computing the orders nu,
    nu+1, ..., nu+n-1. It is also the size of the cyl_i_scal_arr array.
    - z, complex argument of I_nu(z)*exp(-abs(real(z))).
    - cyl_i_scal_arr, array of size n to output I_nu(z)*exp(-abs(real(z))) for
    the orders nu, nu+1, ..., nu+n-1.
    
    Implementation: Similar to the cyl_i_seq() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
void cyl_i_scal_seq(double nu, int n,
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_i_scal_arr) {
    
    cyl_i_full_seq_impl_(nu, n, z, cyl_i_scal_arr, 1);
}

#endif /* BESSEL_LIBRARY_CYL_I_H */