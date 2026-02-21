/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Function pow_dd of the f2c library.
*/

#ifndef BESSEL_LIBRARY_F2C_POW_DD_IMPL_H
#define BESSEL_LIBRARY_F2C_POW_DD_IMPL_H

#include <math.h> /* For math operations and constants */

static inline double f2c_pow_dd_impl_(double *x, double *y) {
  return pow(*x, *y);
}

#endif /* BESSEL_LIBRARY_F2C_POW_DD_IMPL_H */