/*
  Bessel Library: A C library with routines for computing Bessel functions

  File: include/bessel-library/impl/fmax_impl_.h
  Language standards: C99
  References: include/bessel-library/references.txt
  License: include/bessel-library/license.txt
  Repository: <https://github.com/jodesarro/bessel-library>

  Description: Returns the maximum between two real numbers of double type.
*/

#ifndef BESSEL_LIBRARY_FMAX_IMPL_H
#define BESSEL_LIBRARY_FMAX_IMPL_H

static inline double fmax_impl_(double x, double y) { return (x > y ? x : y); }

#endif /* BESSEL_LIBRARY_FMAX_IMPL_H */