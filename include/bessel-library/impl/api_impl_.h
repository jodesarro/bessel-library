/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Defines the API with macros for C++, and for compilation, and for
  header-only or compiled library usage.
*/

#ifndef BESSEL_LIBRARY_API_IMPL_H
#define BESSEL_LIBRARY_API_IMPL_H

#ifdef BESSEL_LIBRARY_EXPORTS_IMPL_
#if defined(_WIN32) || defined(_WIN64)
#define BESSEL_LIBRARY_VISIBILITY_IMPL_ __declspec(dllexport)
#else
#define BESSEL_LIBRARY_VISIBILITY_IMPL_ __attribute__((visibility("default")))
#endif
#elif defined(BESSEL_LIBRARY_IMPORTS)
#if defined(_WIN32) || defined(_WIN64)
#define BESSEL_LIBRARY_VISIBILITY_IMPL_ __declspec(dllimport)
#else
#define BESSEL_LIBRARY_VISIBILITY_IMPL_
#endif
#else
#define BESSEL_LIBRARY_VISIBILITY_IMPL_ static inline
#endif

#ifdef __cplusplus
#define BESSEL_LIBRARY_CPP_IMPL_ extern "C"
#else
#define BESSEL_LIBRARY_CPP_IMPL_
#endif

#define BESSEL_LIBRARY_API_IMPL_                                               \
  BESSEL_LIBRARY_VISIBILITY_IMPL_ BESSEL_LIBRARY_CPP_IMPL_

#endif /* BESSEL_LIBRARY_API_IMPL_H */