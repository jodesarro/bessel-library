/*
  Bessel Library: A C library with routines for computing Bessel functions

  Repository: <https://github.com/jodesarro/bessel-library>
  License: Refer to the LICENSE file in the Repository
  References: Refer to the README.md file in the Repository
  Language standard: C99

  Description: Implements and computes a n-sequency array in double complex type
  for C, or in std::complex<double> for C++, of modified cylindrical Bessel
  functions of the second kind, real order nu, and complex argument z, i.e.,
  {K_nu(z), K_(nu+1)(z), ..., K_(nu+n-1)(z)}, for also negative orders, by means
  of the routines from the Slatec library and recurrence relations for negative
  orders.
*/

#ifndef BESSEL_LIBRARY_CYL_K_FULL_SEQ_IMPL_H
#define BESSEL_LIBRARY_CYL_K_FULL_SEQ_IMPL_H

#include "cplx_c_cpp_impl_.h"
#include "slatec_flags_impl_.h"
#include "slatec_zbesk_impl_.h"
#include <float.h>  /*For DBL_EPSILON */
#include <math.h>   /* For math operations and constants */
#include <stdlib.h> /* For malloc() and free() */

/* Fallback for M_PI if not defined by <math.h> */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
  Implements and computes a n-sequency array of modified cylindrical Bessel
  functions of the second kind, real order nu, and complex argument z, i.e.,
  {K_nu(z), K_(nu+1)(z), ..., K_(nu+n-1)(z)}, for also negative orders, by means
  of the routines from the Slatec library and recurrence relations for negative
  orders.

  Parameters:
  - nu, real order of K_nu(z).
  - n, number n of elements in the sequence for computing the orders nu, nu+1,
  ..., nu+n-1. It is also the size of the cyl_k_arr array.
  - z, complex argument of K_nu(z).
  - cyl_k_arr, array of size n to output K_nu(z) for the orders nu, nu+1, ...,
  nu+n-1
  - scaled, returns the scaled version K_nu(z)*exp(z) if 1.

  Implementation: In general, the implementation is based on the D. E. Amos
  Fortran 77 routines from the Slatec library [3]. Such Fortran routines, and
  all their dependencies, were carefully translated to C. Negative orders are
  handled by Eqs. (6.5.5) of Ref. [2]. When abs(z)=0, it yields INFINITY if
  nu=0, or INFINITY + I * INFINITY otherwise.
*/
static inline void cyl_k_full_seq_impl_(double nu, int n, dcomplex z,
                                        dcomplex *cyl_k_arr, int scaled) {

  int kode = (scaled == 1 ? 2 : 1);
  int nz, ierr;
  double x = creal(z);
  double y = cimag(z);
  double fnu = fabs(nu);
  double nu_m = nu + (double)(n - 1);

  if (cabs(z) < DBL_EPSILON) {

    /* Complex infinity for all orders */
    for (int i = 0; i < n; i++) {
      cyl_k_arr[i] = INFINITY + I * INFINITY;
    }
    if (fabs(nu_m - floor(nu_m)) < DBL_EPSILON && n > fnu) {

      /* Infinity for integer order 0 */
      cyl_k_arr[(int)floor(fnu)] = INFINITY + I * 0.0;
    }

  } else if (nu >= 0.0) {

    /* Positive orders */

    /* Dynamic mem alloc of auxiliary arrays */
    double *ckr = (double *)malloc((n + 1) * sizeof(double));
    double *cki = (double *)malloc((n + 1) * sizeof(double));
    ++ckr;
    ++cki;

    /* Compute zbesk */
    slatec_zbesk_impl_(&x, &y, &fnu, &kode, &n, &ckr[0], &cki[0], &nz, &ierr);
    slatec_flags_zbesk_impl_(ierr, nz);

    /* Store in the array */
    for (int i = 0; i < n; i++) {
      cyl_k_arr[i] = ckr[i] + I * cki[i];
    }

    /* Free auxiliary pointers */
    free(--ckr);
    free(--cki);

  } else if (nu_m <= 0.0) {

    /* Only negative orders */

    double *ckr = (double *)malloc((n + 1) * sizeof(double));
    double *cki = (double *)malloc((n + 1) * sizeof(double));
    ++ckr;
    ++cki;

    /* Compute zbesk */
    double fnu_m = fabs(nu_m);
    slatec_zbesk_impl_(&x, &y, &fnu_m, &kode, &n, &ckr[0], &cki[0], &nz, &ierr);
    slatec_flags_zbesk_impl_(ierr, nz);

    /* Store in the array */
    for (int i = 0; i < n; i++) {
      int tmp = n - 1 - i;
      /* Eq. (6.5.5) of Ref. [2] */
      cyl_k_arr[i] = ckr[tmp] + I * cki[tmp];
    }

    /* Free auxiliary pointers */
    free(--ckr);
    free(--cki);

  } else {

    /* Negative and positive orders */

    int n_m = (int)floor(fabs(nu)) + 1;
    int n_p = n - n_m;
    double fnu_m = fabs(nu + (double)(n_m - 1));

    /* Negative orders */

    /* Dynamic mem alloc of auxiliary arrays */
    double *ckr_m = (double *)malloc((n_m + 1) * sizeof(double));
    double *cki_m = (double *)malloc((n_m + 1) * sizeof(double));
    ++ckr_m;
    ++cki_m;

    /* Compute zbesk */
    slatec_zbesk_impl_(&x, &y, &fnu_m, &kode, &n_m, &ckr_m[0], &cki_m[0], &nz,
                       &ierr);
    slatec_flags_zbesk_impl_(ierr, nz);

    /* Store in the array */
    for (int i = 0; i < n_m; i++) {
      int tmp = n_m - 1 - i;
      /* Eq. (6.5.5) of Ref. [2] */
      cyl_k_arr[i] = ckr_m[tmp] + I * cki_m[tmp];
    }

    /* Free auxiliary pointers */
    free(--ckr_m);
    free(--cki_m);

    /* Positive orders */

    /* Dynamic mem alloc of auxiliary arrays */
    double *ckr_p = (double *)malloc((n_p + 1) * sizeof(double));
    double *cki_p = (double *)malloc((n_p + 1) * sizeof(double));
    ++ckr_p;
    ++cki_p;

    /* Compute zbesk */
    double fnu_p = nu + (double)n_m;
    slatec_zbesk_impl_(&x, &y, &fnu_p, &kode, &n_p, &ckr_p[0], &cki_p[0], &nz,
                       &ierr);
    slatec_flags_zbesk_impl_(ierr, nz);

    /* Store in the array */
    for (int i = n_m; i < n; i++) {
      int tmp1 = i - n_m;
      cyl_k_arr[i] = ckr_p[tmp1] + I * cki_p[tmp1];
    }

    /* Free auxiliary pointers */
    free(--ckr_p);
    free(--cki_p);
  }
}

#endif /* BESSEL_LIBRARY_CYL_K_FULL_SEQ_IMPL_H */