/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/core/cyl_k.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Computes modified cylindrical Bessel functions of the second
        kind, real order, and complex argument.
*/

#ifndef BESSEL_LIBRARY_CYL_K_H
#define BESSEL_LIBRARY_CYL_K_H

#ifndef BESSEL_LIBRARY_STATIC_INLINE_IMPL_
#define BESSEL_LIBRARY_STATIC_INLINE_IMPL_ static inline
#endif

#include "../impl/cplx_c_cpp_impl_.h"
#include "../impl/cyl_k_full_seq_impl_.h"

/*
    Returns the modified cylindrical Bessel function of the second kind, real
    order nu, and complex argument z, i.e., K_nu(z).

    Parameters:
    - nu, real order of K_nu(z).
    - z, complex argument of K_nu(z).
        
    Implementation: Similar to the cyl_k_seq() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdfcplx_impl_ cyl_k(double nu, tpdfcplx_impl_ z) {
    
    /* Array of one size */
    tpdfcplx_impl_ ck[1];
    
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
            
    Implementation: Similar to the cyl_k_seq() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdfcplx_impl_ cyl_k_scal(double nu, tpdfcplx_impl_ z) {
    
    /* Array of one size */
    tpdfcplx_impl_ ck[1];
    
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
    - n, number n of elements in the sequence for computing the orders nu,
    nu+1, ..., nu+n-1. It is also the size of the cyl_k_arr array.
    - z, complex argument of K_nu(z).
    - cyl_k_arr, array of size n to output K_nu(z) for the orders nu,
    nu+1, ..., nu+n-1.
            
    Implementation: In general, the implementation is based on the D. E. Amos
    Fortran 77 routines from the Slatec library [3]. Such Fortran routines,
    and all their dependencies, were carefully translated to C. Negative
    orders are handled by Eqs. (6.5.5) of Ref. [2]. When
    abs(z)=0, it yields INFINITY if nu=0, or INFINITY + I * INFINITY
    otherwise.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
void cyl_k_seq(double nu, int n, tpdfcplx_impl_ z,
    tpdfcplx_impl_ *cyl_k_arr) {

    cyl_k_full_seq_impl_(nu, n, z, cyl_k_arr, 0);
}

/*
    Computes a n-sequency array of scaled versions of modified
    cylindrical Bessel functions of the second kind, real order nu, and
    complex argument z, i.e., {K_nu(z)*exp(z), Y_(nu+1)(z)*exp(z),
    ..., Y_(nu+n-1)(z)*exp(z)}.

    Parameters:
    - nu, real order of K_nu(z)*exp(z).
    - n, number n of elements in the sequence for computing the orders nu,
    nu+1, ..., nu+n-1. It is also the size of the cyl_k_scaled_arr array.
    - z, complex argument of K_nu(z)*exp(z).
    - cyl_k_scaled_arr, array of size n to output K_nu(z)*exp(z) for the
    orders nu, nu+1, ..., nu+n-1.
            
    Implementation: Similar to the cyl_k_seq() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
void cyl_k_scal_seq(double nu, int n,
    tpdfcplx_impl_ z, tpdfcplx_impl_ *cyl_k_scaled_arr) {
    
    cyl_k_full_seq_impl_(nu, n, z, cyl_k_scaled_arr, 1);
}

#endif /* BESSEL_LIBRARY_CYL_K_H */