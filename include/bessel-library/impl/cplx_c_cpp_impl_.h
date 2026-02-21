/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99 + C++98

  Description: Defines macros and typedefs for ensuring C++ compatibility.
*/

#ifndef BESSEL_LIBRARY_CPLX_C_CPP_IMPL_H
#define BESSEL_LIBRARY_CPLX_C_CPP_IMPL_H

#ifdef __cplusplus

/* Includes, typedefs and/or macros for C++98 compatibility */

#include <complex> /* For complex numbers */
#define I_IMPL_ std::complex<double>(0.0, 1.0)
#define CPLX_IMPL_(x, y) std::complex<double>(x, y)
#define creal(z) std::real(z)
#define cimag(z) std::imag(z)
#define cabs(z) std::abs(z)
#define cexp(z) std::exp(z)
#define sin(z) std::sin(z)
#define cos(z) std::cos(z)
typedef std::complex<double> tpdfcplx_impl_;

#else

/* C99 */

#include <complex.h> /* For complex numbers */
#define I_IMPL_ I
#define CPLX_IMPL_(x, y) (x + I * y)
typedef double complex tpdfcplx_impl_;

#endif /* __cplusplus */

#endif /* BESSEL_LIBRARY_CPLX_C_CPP_IMPL_H */