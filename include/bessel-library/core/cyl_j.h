/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Computes cylindrical Bessel functions of the first kind, real
  order, and complex argument.
*/

#ifndef BESSEL_LIBRARY_CYL_J_H
#define BESSEL_LIBRARY_CYL_J_H

#include "../impl/api_impl_.h"
#include "../impl/cplx_c_cpp_impl_.h"

#ifndef BESSEL_LIBRARY_IMPORTS
#include "../impl/cyl_j_full_seq_impl_.h"
#endif

/*
  Returns the cylindrical Bessel function of the first kind, real order nu, and
  complex argument z, i.e., J_nu(z).

  Parameters:
  - nu, real order of J_nu(z).
  - z, complex argument of J_nu(z).

  Implementation: Similar to the cyl_j_seq() function.
*/
BESSEL_LIBRARY_API_IMPL_
dcomplex cyl_j(double nu, dcomplex z)
#ifndef BESSEL_LIBRARY_IMPORTS
{
  /* Array of one size */
  dcomplex cj[1];

  /* Compute cyl_j_full_seq_impl_ */
  cyl_j_full_seq_impl_(nu, 1, z, cj, 0);

  /* Return */
  return cj[0];
}
#else
    ;
#endif

/*
  Returns the scaled version of the cylindrical Bessel function of the first
  kind, real order nu, and complex argument z, i.e., J_nu(z)*exp(-abs(imag(z))).

  Parameters:
  - nu, real order of J_nu(z)*exp(-abs(imag(z))).
  - z, complex argument of J_nu(z)*exp(-abs(imag(z))).

  Implementation: Similar to the cyl_j_seq() function.
*/
BESSEL_LIBRARY_API_IMPL_
dcomplex cyl_j_scal(double nu, dcomplex z)
#ifndef BESSEL_LIBRARY_IMPORTS
{
  /* Array of one size */
  dcomplex cj[1];

  /* Compute cyl_j_full_seq_impl_ */
  cyl_j_full_seq_impl_(nu, 1, z, cj, 1);

  /* Return */
  return cj[0];
}
#else
    ;
#endif

/*
  Computes a n-sequency array of cylindrical Bessel functions of the first kind,
  real order nu, and complex argument z, i.e., {J_nu(z), J_(nu+1)(z), ...,
  J_(nu+n-1)(z)}.

  Parameters:
  - nu, real order of J_nu(z).
  - n, number n of elements in the sequence for computing the orders nu, nu+1,
  ..., nu+n-1. It is also the size of the cyl_j_arr array.
  - z, complex argument of J_nu(z).
  - cyl_j_arr, array of size n to output J_nu(z) for the orders nu, nu+1, ...,
  nu+n-1

  Implementation: In general, the implementation is based on the D. E. Amos
  Fortran 77 routines from the Slatec library [3]. Such Fortran routines, and
  all their dependencies, were carefully translated to C. Negative orders are
  handled by Eqs. (5.4.2) and (5.5.4) of Ref. [2] for, respectively, nu integer
  and nu real; in the latter case, it yields INFINITY + I * INFINITY when
  abs(z)=0.
*/
BESSEL_LIBRARY_API_IMPL_
void cyl_j_seq(double nu, int n, dcomplex z, dcomplex *cyl_j_arr)
#ifndef BESSEL_LIBRARY_IMPORTS
{
  cyl_j_full_seq_impl_(nu, n, z, cyl_j_arr, 0);
}
#else
    ;
#endif

/*
  Computes a n-sequency array of scaled versions of cylindrical Bessel functions
  of the first kind, real order nu, and complex argument z, i.e.,
  {J_nu(z)*exp(-abs(imag(z))), J_(nu+1)(z)*exp(-abs(imag(z))), ...,
  J_(nu+n-1)(z)*exp(-abs(imag(z)))}.

  Parameters:
  - nu, real order of J_nu(z)*exp(-abs(imag(z))).
  - n, number n of elements in the sequence for computing the orders nu, nu+1,
  ..., nu+n-1. It is also the size of the cyl_j_scal_arr array.
  - z, complex argument of J_nu(z)*exp(-abs(imag(z))).
  - cyl_j_scal_arr, array of size n to output J_nu(z)*exp(-abs(imag(z))) for the
  orders nu, nu+1, ..., nu+n-1

  Implementation: Similar to the cyl_j_seq() function.
*/
BESSEL_LIBRARY_API_IMPL_
void cyl_j_scal_seq(double nu, int n, dcomplex z, dcomplex *cyl_j_scal_arr)
#ifndef BESSEL_LIBRARY_IMPORTS
{
  cyl_j_full_seq_impl_(nu, n, z, cyl_j_scal_arr, 1);
}
#else
    ;
#endif

#endif /* BESSEL_LIBRARY_CYL_J_H */