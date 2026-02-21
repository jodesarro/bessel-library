/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Returns the Airy function of the first kind and complex argument
  z, i.e., Ai(z), in double complex type for C, or in std::complex<double> type
  for C++, by means of the routines from the Slatec library.
*/

#ifndef BESSEL_LIBRARY_AIRY_AI_IMPL_H
#define BESSEL_LIBRARY_AIRY_AI_IMPL_H

#include "cplx_c_cpp_impl_.h"
#include "slatec_flags_impl_.h"
#include "slatec_zairy_impl_.h"

/*
  Returns the Airy function of the first kind and complex argument z, i.e.,
  Ai(z) by means of the routines from the Slatec library.

  Parameters:
  - z, complex argument of Ai(z).
  - derivative, replaces Ai(z) by dAi(z)/dz if 1.
  - scaled, returns the scaled version Ai(z)*exp((2/3)*pow(z,3/2)) if 1.

  Implementation: In general, the implementation is based on the D. E. Amos
  Fortran 77 routines from the Slatec library [3]. Such Fortran routines, and
  all their dependencies, were carefully translated to C.
*/
static inline tpdfcplx_impl_ airy_ai_impl_(tpdfcplx_impl_ z, int derivative,
                                           int scaled) {

  int kode = (scaled == 1 ? 2 : 1);
  int nz, ierr;
  double x = creal(z);
  double y = cimag(z);

  /* Auxiliary arrays */
  double air_arr[2] = {0.0, 0.0}, aii_arr[2] = {0.0, 0.0};

  /* Compute zairy */
  slatec_zairy_impl_(&x, &y, &derivative, &kode, &air_arr[1], &aii_arr[1], &nz,
                     &ierr);
  slatec_flags_zairy_impl_(ierr, nz);

  /* Return */
  return air_arr[1] + I_IMPL_ * aii_arr[1];
}

#endif /* BESSEL_LIBRARY_AIRY_AI_IMPL_H */