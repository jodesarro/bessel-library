/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/core/airy_bi.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 with guards for C++98 compatibility
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Returns, in double complex type for C, or in std::complex<double> type
        for C++, the Airy functions of the second kind and complex argument.
*/

#ifndef BESSEL_LIBRARY_AIRY_BI_H
#define BESSEL_LIBRARY_AIRY_BI_H

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

#include "../impl/airy_bi_impl_.h"

/*
    Returns the Airy function of the second kind and complex argument z, i.e.,
    Bi(z).

    Parameter:
    - z, complex argument of Bi(z).

    Implementation: In general, the implementation is based on the D. E. Amos
    Fortran 77 routines of the Slatec library [3] Such Fortran routines,
    and all their dependencies, were carefully translated to C.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdcomplex_impl_ airy_bi(tpdcomplex_impl_ z) {
    return airy_bi_impl_(z, 0, 0);
}

/*
    Returns the first derivative of the Airy function of the second kind and
    complex argument z, i.e., dBi(z)/dz.

    Parameter:
    - z, complex argument of dBi(z)/dz.
        
    Implementation: Similar to the airy_bi() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdcomplex_impl_ airy_bi_diff(tpdcomplex_impl_ z) {
    return airy_bi_impl_(z, 1, 0);
}

/*
    Returns the scaled version of the Airy function of the second kind and
    complex argument z, i.e., Bi(z)*exp(-abs(real((2/3)*pow(z,3/2)))).

    Parameter:
    - z, complex argument of Bi(z)*exp(-abs(real((2/3)*pow(z,3/2)))).
        
    Implementation: Similar to the airy_bi() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdcomplex_impl_ airy_bi_scal(tpdcomplex_impl_ z) {
    return airy_bi_impl_(z, 0, 1);
}

/*
    Returns the scaled version of the first derivative of the Airy function of
    the second kind and complex argument z, i.e.,
    (dBi(z)/dz)*exp(-abs(real((2/3)*pow(z,3/2)))).

    Parameter:
    - z, complex argument of (dBi(z)/dz)*exp(-abs(real((2/3)*pow(z,3/2)))).
    
    Implementation: Similar to the airy_bi() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdcomplex_impl_ airy_bi_diff_scal(tpdcomplex_impl_ z) {
    return airy_bi_impl_(z, 1, 1);
}

#endif /* BESSEL_LIBRARY_AIRY_BI_H */