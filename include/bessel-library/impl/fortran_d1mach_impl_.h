/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Fortran's D1MACH function.
*/

#ifndef BESSEL_LIBRARY_FORTRAN_D1MACH_IMPL_H
#define BESSEL_LIBRARY_FORTRAN_D1MACH_IMPL_H

#include <float.h> /* For DBL limits */

static inline double fortran_d1mach_impl_(int *c_un) {
  switch (*c_un) {
  case 1:
    return DBL_MIN;
    break;
  case 2:
    return DBL_MAX;
    break;
  case 4:
    return DBL_EPSILON;
    break;
  case 5:
    return 0.3010299956639811952137; /* log10(2) */
    break;
  default:
    return 0.0;
    break;
  }
}

#endif /* BESSEL_LIBRARY_FORTRAN_D1MACH_IMPL_H */