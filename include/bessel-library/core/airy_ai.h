/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: bessel-library/core/airy_ai.h
    Version: bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: bessel-library/references.txt

    Description:
        Returns the Airy functions of the first kind and complex argument.
*/

#ifndef BESSEL_LIBRARY_AIRY_AI_H
#define BESSEL_LIBRARY_AIRY_AI_H

#ifndef BESSEL_LIBRARY_STATIC_INLINE_IMPL_
#define BESSEL_LIBRARY_STATIC_INLINE_IMPL_ static inline
#endif

#include "../impl/cplx_c_cpp_impl_.h"
#include "../impl/airy_ai_impl_.h"

/*
    Returns the Airy function of the first kind and complex argument z, i.e.,
    Ai(z).

    Parameter:
    - z, complex argument of Ai(z).

    Implementation: In general, the implementation is based on the D. E.
    Amos Fortran 77 routines from the Slatec library [3]. Such Fortran
    routines, and all their dependencies, were carefully translated to C.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdfcplx_impl_ airy_ai(tpdfcplx_impl_ z) {
    return airy_ai_impl_(z, 0, 0);
}

/*
    Returns the first derivative of the Airy function of the first kind and
    complex argument z, i.e., dAi(z)/dz.

    Parameter:
    - z, complex argument of dAi(z)/dz.
    
    Implementation: Similar to the airy_ai() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdfcplx_impl_ airy_ai_diff(tpdfcplx_impl_ z) {
    return airy_ai_impl_(z, 1, 0);
}

/*
    Returns the scaled version of the Airy function of the first kind and
    complex argument z, i.e., Ai(z)*exp((2/3)*pow(z,3/2)).

    Parameter:
    - z, complex argument of Ai(z)*exp((2/3)*pow(z,3/2)).

    Implementation: Similar to the airy_ai() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdfcplx_impl_ airy_ai_scal(tpdfcplx_impl_ z) {
    return airy_ai_impl_(z, 0, 1);
}

/*
    Returns the scaled version of the first derivative of the Airy function of
    the first kind and complex argument z, i.e.,
    (dAi(z)/dz)*exp((2/3)*pow(z,3/2)).

    Parameter:
    - z, complex argument of (dAi(z)/dz)*exp((2/3)*pow(z,3/2)).

    Implementation: Similar to the airy_ai() function.
*/
BESSEL_LIBRARY_STATIC_INLINE_IMPL_
tpdfcplx_impl_ airy_ai_diff_scal(tpdfcplx_impl_ z) {
    return airy_ai_impl_(z, 1, 1);
}

#endif /* BESSEL_LIBRARY_AIRY_AI_H */