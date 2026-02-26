/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99 and C++98

  Description: Defines macros and typedefs for ensuring C++ compatibility.
*/

#ifndef BESSEL_LIBRARY_CPLX_C_CPP_IMPL_H
#define BESSEL_LIBRARY_CPLX_C_CPP_IMPL_H

#ifdef __cplusplus

/* Includes, typedefs and/or macros for C++98 compatibility */

#include <complex>                     /* For complex numbers */
typedef std::complex<double> dcomplex; /* std::complex<double> */
#define I (dcomplex(0.0, 1.0))         /* std::complex<double>(0.0, 1.0) */
#define creal(z) (std::real(z))        /* std::real(z) */
#define cimag(z) (std::imag(z))        /* std::imag(z) */
#define cabs(z) (std::abs(z))          /* std::abs(z) */
#define cexp(z) (std::exp(z))          /* std::exp(z) */

#else

/* C99 */

#include <complex.h> /* For complex numbers */
typedef double complex dcomplex; /* double complex */

#endif /* __cplusplus */

#endif /* BESSEL_LIBRARY_CPLX_C_CPP_IMPL_H */