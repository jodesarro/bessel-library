/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/core/cyl_h1.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Computes cylindrical Hankel functions of the first kind, real
        order, and complex argument.
*/

#ifndef BESSEL_LIBRARY_CYL_H1_H
#define BESSEL_LIBRARY_CYL_H1_H

#include "../impl/api_impl_.h"
#include "../impl/cplx_c_cpp_impl_.h"

#ifndef BESSEL_LIBRARY_IMPORTS
#include "../impl/cyl_h1_full_seq_impl_.h"
#endif

/*
    Returns the cylindrical Hankel function of the first kind, real order nu,
    and complex argument z, i.e., H1_nu(z).

    Parameters:
    - nu, real order of H1_nu(z).
    - z, complex argument of H1_nu(z).
            
    Implementation: Similar to the cyl_h1_seq() function.
*/
BESSEL_LIBRARY_API_IMPL_
tpdfcplx_impl_ cyl_h1(double nu, tpdfcplx_impl_ z)
#ifndef BESSEL_LIBRARY_IMPORTS
{
    /* Array of one size */
    tpdfcplx_impl_ ch1[1];
    
    /* Compute cyl_h1_full_seq_impl_ */
    cyl_h1_full_seq_impl_(nu, 1, z, ch1, 0);

    /* Return */
    return ch1[0];
}
#else
;
#endif

/*
    Returns the scaled version of the cylindrical Hankel function of the
    first kind, real order nu, and complex argument z, i.e.,
    H1_nu(z)*exp(-i*z).

    Parameters:
    - nu, real order of H1_nu(z)*exp(-i*z).
    - z, complex argument of H1_nu(z)*exp(-i*z).
                
    Implementation: Similar to the cyl_h1_seq() function.
*/
BESSEL_LIBRARY_API_IMPL_
tpdfcplx_impl_ cyl_h1_scal(double nu, tpdfcplx_impl_ z)
#ifndef BESSEL_LIBRARY_IMPORTS
{
    /* Array of one size */
    tpdfcplx_impl_ ch1[1];
    
    /* Compute cyl_h1_full_seq_impl_ */
    cyl_h1_full_seq_impl_(nu, 1, z, ch1, 1);

    /* Return */
    return ch1[0];
}
#else
;
#endif

/*
    Computes a n-sequency array of cylindrical Hankel functions of the first
    kind, real order nu, and complex argument z, i.e., {H1_nu(z),
    H1_(nu+1)(z), ..., H1_(nu+n-1)(z)}.

    Parameters:
    - nu, real order of H1_nu(z).
    - n, number n of elements in the sequence for computing the orders nu,
    nu+1, ..., nu+n-1. It is also the size of the cyl_h1_arr array.
    - z, complex argument of H1_nu(z).
    - cyl_h1_arr, array of size n to output H1_nu(z) for the orders nu, nu+1,
    ..., nu+n-1.

    Implementation: In general, the implementation is based on the D. E. Amos
    Fortran 77 routines from the Slatec library [3]. Such Fortran routines,
    and all their dependencies, were carefully translated to C. Negative
    orders are handled by Eq. (9.1.6) of Ref. [1].
    It yields INFINITY + I * INFINITY when abs(z)=0.
*/
BESSEL_LIBRARY_API_IMPL_
void cyl_h1_seq(double nu, int n,
    tpdfcplx_impl_ z, tpdfcplx_impl_ *cyl_h1_arr)
#ifndef BESSEL_LIBRARY_IMPORTS
{
    cyl_h1_full_seq_impl_(nu, n, z, cyl_h1_arr, 0);
}
#else
;
#endif

/*
    Computes a n-sequency array of scaled versions of cylindrical
    Hankel functions of the first kind, real order nu, and complex argument z,
    i.e., {H1_nu(z)*exp(-i*z), H1_(nu+1)(z)*exp(-i*z), ...,
    H1_(nu+n-1)(z)*exp(-i*z)}.

    Parameters:
    - nu, real order of H1_nu(z)*exp(-i*z).
    - n, number n of elements in the sequence for computing the orders nu,
    nu+1, ..., nu+n-1. It is also the size of the cyl_h2_scal_arr array.
    - z, complex argument of H1_nu(z)*exp(-i*z).
    - cyl_h2_scal_arr, array of size n to output  H1_nu(z)*exp(-i*z) for the
    orders nu, nu+1, ..., nu+n-1.
                
    Implementation: Similar to the cyl_h1_seq() function.
*/
BESSEL_LIBRARY_API_IMPL_
void cyl_h1_scal_seq(double nu, int n,
    tpdfcplx_impl_ z, tpdfcplx_impl_ *cyl_h2_scal_arr)
#ifndef BESSEL_LIBRARY_IMPORTS
{
    cyl_h1_full_seq_impl_(nu, n, z, cyl_h2_scal_arr, 1);
}
#else
;
#endif

#endif /* BESSEL_LIBRARY_CYL_H1_H */