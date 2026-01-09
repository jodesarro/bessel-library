/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: bessel-library/core/airy_ai.h
    Version: bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 with guards for C++98 compatibility
    References: bessel-library/references.txt

    Description:
        Returns, in double complex type for C, or in std::complex<double> type
        for C++, the Airy functions of the first kind and complex argument.
*/

#ifndef BESSEL_LIBRARY_AIRY_AI_H
#define BESSEL_LIBRARY_AIRY_AI_H

#ifdef __cplusplus

/* Includes, typedefs and/or macros for C++98 compatibility */
#include <complex> /* for complex numbers */
typedef std::complex<double> tpdcomplex_impl_;

#else

#include <complex.h> /* for complex numbers */
typedef double complex tpdcomplex_impl_;

#endif /* __cplusplus */

#include "../impl/airy_ai_impl_.h"

/*
    Returns the Airy function of the first kind and complex argument z, i.e.,
    Ai(z).

    Parameter:
    - z, complex argument of Ai(z).
*/
static inline tpdcomplex_impl_ airy_ai(tpdcomplex_impl_ z) {
    return airy_ai_impl_(z, 0, 0);
}

/*
    Returns the first derivative of the Airy function of the first kind and
    complex argument z, i.e., dAi(z)/dz.

    Parameter:
    - z, complex argument of dAi(z)/dz.
*/
static inline tpdcomplex_impl_ airy_ai_deriv(tpdcomplex_impl_ z) {
    return airy_ai_impl_(z, 1, 0);
}

/*
    Returns the scaled version of the Airy function of the first kind and
    complex argument z, i.e., Ai(z)*exp((2/3)*pow(z,3/2)).

    Parameter:
    - z, complex argument of Ai(z)*exp((2/3)*pow(z,3/2)).
*/
static inline tpdcomplex_impl_ airy_ai_scal(tpdcomplex_impl_ z) {
    return airy_ai_impl_(z, 0, 1);
}

/*
    Returns the scaled version of the first derivative of the Airy function of
    the first kind and complex argument z, i.e.,
    (dAi(z)/dz)*exp((2/3)*pow(z,3/2)).

    Parameter:
    - z, complex argument of (dAi(z)/dz)*exp((2/3)*pow(z,3/2)).
*/
static inline tpdcomplex_impl_ airy_ai_deriv_scal(tpdcomplex_impl_ z) {
    return airy_ai_impl_(z, 1, 1);
}

#endif /* BESSEL_LIBRARY_AIRY_AI_H */