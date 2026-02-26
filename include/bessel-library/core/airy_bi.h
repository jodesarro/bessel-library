/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description:Returns the Airy functions of the second kind and complex
  argument.
*/

#ifndef BESSEL_LIBRARY_AIRY_BI_H
#define BESSEL_LIBRARY_AIRY_BI_H

#include "../impl/api_impl_.h"
#include "../impl/cplx_c_cpp_impl_.h"

#ifndef BESSEL_LIBRARY_IMPORTS
#include "../impl/airy_bi_impl_.h"
#endif

/*
  Returns the Airy function of the second kind and complex argument z, i.e.,
  Bi(z).

  Parameter:
  - z, complex argument of Bi(z).

  Implementation: In general, the implementation is based on the D. E. Amos
  Fortran 77 routines from the Slatec library [3] Such Fortran routines, and all
  their dependencies, were carefully translated to C.
*/
BESSEL_LIBRARY_API_IMPL_
dcomplex airy_bi(dcomplex z)
#ifndef BESSEL_LIBRARY_IMPORTS
{
  return airy_bi_impl_(z, 0, 0);
}
#else
    ;
#endif

/*
  Returns the first derivative of the Airy function of the second kind and
  complex argument z, i.e., dBi(z)/dz.

  Parameter:
  - z, complex argument of dBi(z)/dz.

  Implementation: Similar to the airy_bi() function.
*/
BESSEL_LIBRARY_API_IMPL_
dcomplex airy_bi_diff(dcomplex z)
#ifndef BESSEL_LIBRARY_IMPORTS
{
  return airy_bi_impl_(z, 1, 0);
}
#else
    ;
#endif

/*
  Returns the scaled version of the Airy function of the second kind and complex
  argument z, i.e., Bi(z)*exp(-abs(real((2/3)*pow(z,3/2)))).

  Parameter:
  - z, complex argument of Bi(z)*exp(-abs(real((2/3)*pow(z,3/2)))).

  Implementation: Similar to the airy_bi() function.
*/
BESSEL_LIBRARY_API_IMPL_
dcomplex airy_bi_scal(dcomplex z)
#ifndef BESSEL_LIBRARY_IMPORTS
{
  return airy_bi_impl_(z, 0, 1);
}
#else
    ;
#endif

/*
  Returns the scaled version of the first derivative of the Airy function of the
  second kind and complex argument z, i.e.,
  (dBi(z)/dz)*exp(-abs(real((2/3)*pow(z,3/2)))).

  Parameter:
  - z, complex argument of (dBi(z)/dz)*exp(-abs(real((2/3)*pow(z,3/2)))).

  Implementation: Similar to the airy_bi() function.
*/
BESSEL_LIBRARY_API_IMPL_
dcomplex airy_bi_diff_scal(dcomplex z)
#ifndef BESSEL_LIBRARY_IMPORTS
{
  return airy_bi_impl_(z, 1, 1);
}
#else
    ;
#endif

#endif /* BESSEL_LIBRARY_AIRY_BI_H */