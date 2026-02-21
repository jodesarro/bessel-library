/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Returns the minimum between two integer numbers of int type.
*/

#ifndef BESSEL_LIBRARY_MIN_IMPL_H
#define BESSEL_LIBRARY_MIN_IMPL_H

static inline int min_impl_(int x, int y) { return (x < y ? x : y); }

#endif /* BESSEL_LIBRARY_MIN_IMPL_H */