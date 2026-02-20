/*
  Bessel Library: A C library with routines for computing Bessel functions

  File: include/bessel-library/impl/max_impl_.h
  Language standards: C99
  References: include/bessel-library/references.txt
  License: include/bessel-library/license.txt
  Repository: <https://github.com/jodesarro/bessel-library>

  Description: Returns the maximum between two integer numbers of int type.
*/

#ifndef BESSEL_LIBRARY_MAX_IMPL_H
#define BESSEL_LIBRARY_MAX_IMPL_H

static inline int max_impl_(int x, int y) { return (x > y ? x : y); }

#endif /* BESSEL_LIBRARY_MAX_IMPL_H */