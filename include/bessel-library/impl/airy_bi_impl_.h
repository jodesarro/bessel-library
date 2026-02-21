/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Returns the Airy function of the second kind and complex argument
  z, i.e., Bi(z), in double complex type for C, or in std::complex<double> type
  for C++, by means of the routines from the Slatec library.
*/

#ifndef BESSEL_LIBRARY_AIRY_BI_IMPL_H
#define BESSEL_LIBRARY_AIRY_BI_IMPL_H

#include "cplx_c_cpp_impl_.h"
#include "slatec_flags_impl_.h"
#include "slatec_zbiry_impl_.h"

/*
  Returns the Airy function of the second kind and complex argument z, i.e.,
  Bi(z), by means of the routines from the Slatec library.

  Parameters:
  - z, complex argument of Bi(z).
  - derivative, replaces Bi(z) by dBi(z)/dz if 1.
  - scaled, returns the scaled version Bi(z)*exp(abs(real((2/3)*pow(z,3/2))))
  if 1.

  Implementation: In general, the implementation is based on the D. E. Amos
  Fortran 77 routines from the Slatec library [3] Such Fortran routines, and all
  their dependencies, were carefully translated to C.
*/
static inline tpdfcplx_impl_ airy_bi_impl_(tpdfcplx_impl_ z, int derivative,
                                           int scaled) {

  int kode = (scaled == 1 ? 2 : 1);
  int ierr;
  double x = creal(z);
  double y = cimag(z);

  /* Auxiliary arrays */
  double bir_arr[2] = {0.0, 0.0}, bii_arr[2] = {0.0, 0.0};

  /* Compute zbiry */
  slatec_zbiry_impl_(&x, &y, &derivative, &kode, &bir_arr[1], &bii_arr[1],
                     &ierr);
  slatec_flags_zbiry_impl_(ierr);

  /* Return */
  return bir_arr[1] + I_IMPL_ * bii_arr[1];
}

#endif /* BESSEL_LIBRARY_AIRY_BI_IMPL_H */