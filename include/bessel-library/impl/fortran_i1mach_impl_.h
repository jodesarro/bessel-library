/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Fortran's I1MACH function.
*/

#ifndef BESSEL_LIBRARY_FORTRAN_I1MACH_IMPL_H
#define BESSEL_LIBRARY_FORTRAN_I1MACH_IMPL_H

#include <float.h>  /* For DBL limits */
#include <limits.h> /* For INT_MAX */

static inline int fortran_i1mach_impl_(int *c_un) {
  switch (*c_un) {
  case 9:
    return INT_MAX;
    break;
  case 14:
    return DBL_MANT_DIG;
    break;
  case 15:
    return -DBL_MIN_EXP;
    break;
  case 16:
    return DBL_MAX_EXP;
    break;
  default:
    return 0;
    break;
  }
}

#endif /* BESSEL_LIBRARY_FORTRAN_I1MACH_IMPL_H */