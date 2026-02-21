/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Returns the minimum between two real numbers of double type.
*/

#ifndef BESSEL_LIBRARY_FMIN_IMPL_H
#define BESSEL_LIBRARY_FMIN_IMPL_H

static inline double fmin_impl_(double x, double y) { return (x < y ? x : y); }

#endif /* BESSEL_LIBRARY_FMIN_IMPL_H */