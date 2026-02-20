/*
  Bessel Library: A C library with routines for computing Bessel functions

  File: include/bessel-library/impl/min_impl.h
  Language standards: C99
  References: include/bessel-library/references.txt
  License: include/bessel-library/license.txt
  Repository: <https://github.com/jodesarro/bessel-library>

  Description: Returns the minimum between two integer numbers of int type.
*/

#ifndef BESSEL_LIBRARY_MIN_IMPL_H
#define BESSEL_LIBRARY_MIN_IMPL_H

static inline int min_impl_(int x, int y) { return (x < y ? x : y); }

#endif /* BESSEL_LIBRARY_MIN_IMPL_H */