/*
  Bessel Library: A C library with routines for computing Bessel functions

  File: include/bessel-library/impl/f2c_d_sign_impl_.h
  Language standards: C99
  References: include/bessel-library/references.txt
  License: include/bessel-library/license.txt
  Repository: <https://github.com/jodesarro/bessel-library>

  Description: Function d_sign of the f2c library.
*/

#ifndef BESSEL_LIBRARY_F2C_D_SIGN_IMPL_H
#define BESSEL_LIBRARY_F2C_D_SIGN_IMPL_H

#include <math.h> /* For math operations and constants */

static inline double f2c_d_sign_impl_(double *x, double *y) {
  return (*y >= 0. ? fabs(*x) : -fabs(*x));
}

#endif /* BESSEL_LIBRARY_F2C_D_SIGN_IMPL_H */