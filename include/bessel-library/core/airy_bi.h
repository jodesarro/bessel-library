/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/core/airy_bi.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Returns the Airy functions of the second kind and complex argument.
*/

#ifndef BESSEL_LIBRARY_AIRY_BI_H
#define BESSEL_LIBRARY_AIRY_BI_H

#ifndef BESSEL_LIBRARY_STATIC_INLINE_IMPL_
#define BESSEL_LIBRARY_STATIC_INLINE_IMPL_ static inline
#endif

#include "../impl/cplx_c_cpp_impl_.h"
#include "../impl/airy_bi_impl_.h"

/*
    Returns the Airy function of the second kind and complex argument z, i.e.,
    Bi(z).

    Parameter:
    - z, complex argument of Bi(z).

    Implementation: In general, the implementation is based on the D. E. Amos
    Fortran 77 routines from the Slatec library [3] Such Fortran routines,
    and all their dependencies, were carefully translated to C.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdfcplx_impl_ airy_bi(tpdfcplx_impl_ z) {
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
tpdfcplx_impl_ airy_bi_diff(tpdfcplx_impl_ z) {
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
tpdfcplx_impl_ airy_bi_scal(tpdfcplx_impl_ z) {
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
tpdfcplx_impl_ airy_bi_diff_scal(tpdfcplx_impl_ z) {
    return airy_bi_impl_(z, 1, 1);
}

#endif /* BESSEL_LIBRARY_AIRY_BI_H */