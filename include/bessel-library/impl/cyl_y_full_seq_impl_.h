/*
  Bessel Library: A C library with routines for computing Bessel functions

  File: include/bessel-library/impl/cyl_y_full_seq_impl_.h
  Language standards: C99
  References: include/bessel-library/references.txt
  License: include/bessel-library/license.txt
  Repository: <https://github.com/jodesarro/bessel-library>

  Description: Implements and computes a n-sequency array in double complex type
  for C, or std::complex<double> for C++, of cylindrical Bessel functions of the
  second kind, real order nu, and complex argument z, i.e., {Y_nu(z),
  Y_(nu+1)(z), ..., Y_(nu+n-1)(z)}, for also negative orders, by means of the
  routines from the Slatec library and recurrence relations for negative orders.
*/

#ifndef BESSEL_LIBRARY_CYL_Y_FULL_SEQ_IMPL_H
#define BESSEL_LIBRARY_CYL_Y_FULL_SEQ_IMPL_H

#include "cplx_c_cpp_impl_.h"
#include "slatec_flags_impl_.h"
#include "slatec_zbesj_impl_.h"
#include "slatec_zbesy_impl_.h"
#include <float.h>  /* For DBL_EPSILON */
#include <math.h>   /* For math operations and constants */
#include <stdlib.h> /* For malloc() and free() */

/* Fallback for M_PI if not defined by <math.h> */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
  Implements and computes a n-sequency array of cylindrical Bessel functions of
  the second kind, real order nu, and complex argument z, i.e., {Y_nu(z),
  Y_(nu+1)(z), ..., Y_(nu+n-1)(z)}, for also negative orders, by means of the
  routines from the Slatec library and recurrence relations for negative orders.

  Parameters:
  - nu, real order of Y_nu(z).
  - n, number n of elements in the sequence for computing the orders nu, nu+1,
  ..., nu+n-1. It is also the size of the cyl_y_arr array.
  - z, complex argument of Y_nu(z).
  - cyl_y_arr, array of size n to output Y_nu(z) for the orders nu, nu+1, ...,
  nu+n-1.
  - scaled, returns the scaled version Y_nu(z)*exp(-abs(imag(z))) if 1.

  Implementation: In general, the implementation is based on the D. E. Amos
  Fortran 77 routines from the Slatec library [3]. Such Fortran routines, and
  all their dependencies, were carefully translated to C. Negative orders are
  handled by Eqs. (5.4.2) and (5.5.4) of Ref. [2] for, respectively, nu integer
  and nu real. When abs(z)=0, it yields -INFINITY if nu=0, or INFINITY + I *
  INFINITY otherwise.
*/
static inline void cyl_y_full_seq_impl_(double nu, int n, tpdfcplx_impl_ z,
                                        tpdfcplx_impl_ *cyl_y_arr, int scaled) {

  int kode = (scaled == 1 ? 2 : 1);
  int nz, ierr;
  double x = creal(z);
  double y = cimag(z);
  double fnu = fabs(nu);
  double nu_m = nu + (double)(n - 1);

  if (cabs(z) < DBL_EPSILON) {

    /* Store in the array */
    for (int i = 0; i < n; i++) {
      /* Complex infinity for i!=0 orders */
      cyl_y_arr[i] = CPLX_IMPL_(INFINITY, INFINITY);
    }
    if (fabs(nu_m - floor(nu_m)) < DBL_EPSILON && n > fnu) {
      /* Complex infinity for i=0 order */
      cyl_y_arr[(int)floor(fnu)] = CPLX_IMPL_(-INFINITY, 0.0);
    }

  } else if (nu >= 0.0) {

    /* Only positive orders, including the order 0 */

    /* Dynamic mem alloc of auxiliary arrays */
    double *cyr_ptr = (double *)malloc((n + 1) * sizeof(double));
    double *cyi_ptr = (double *)malloc((n + 1) * sizeof(double));
    double *cwrkr_ptr = (double *)malloc((n + 1) * sizeof(double));
    double *cwrki_ptr = (double *)malloc((n + 1) * sizeof(double));
    ++cyr_ptr;
    ++cyi_ptr;
    ++cwrkr_ptr;
    ++cwrki_ptr;

    /* Compute zbesy */
    slatec_zbesy_impl_(&x, &y, &fnu, &kode, &n, &cyr_ptr[0], &cyi_ptr[0], &nz,
                       &cwrkr_ptr[0], &cwrki_ptr[0], &ierr);
    slatec_flags_zbesy_impl_(ierr, nz);

    /* Store in the array */
    for (int i = 0; i < n; i++) {
      cyl_y_arr[i] = cyr_ptr[i] + I_IMPL_ * cyi_ptr[i];
    }

    /* Free auxiliary pointers */
    free(--cyr_ptr);
    free(--cyi_ptr);
    free(--cwrkr_ptr);
    free(--cwrki_ptr);

  } else if (nu_m <= 0.0) {

    /* Only negative orders, including the order 0 */

    if (fabs(floor(fnu) - fnu) < DBL_EPSILON) {

      /* Integer negative orders */

      /* Dynamic mem alloc of auxiliary arrays */
      double *cyr_ptr = (double *)malloc((n + 1) * sizeof(double));
      double *cyi_ptr = (double *)malloc((n + 1) * sizeof(double));
      double *cwrkr_ptr = (double *)malloc((n + 1) * sizeof(double));
      double *cwrki_ptr = (double *)malloc((n + 1) * sizeof(double));
      ++cyr_ptr;
      ++cyi_ptr;
      ++cwrkr_ptr;
      ++cwrki_ptr;

      /* Compute zbesy */
      double fnu_m = fabs(nu_m);
      slatec_zbesy_impl_(&x, &y, &fnu_m, &kode, &n, &cyr_ptr[0], &cyi_ptr[0],
                         &nz, &cwrkr_ptr[0], &cwrki_ptr[0], &ierr);
      slatec_flags_zbesy_impl_(ierr, nz);

      /* Store in the array */
      for (int i = 0; i < n; i++) {
        int tmp = n - 1 - i;
        cyl_y_arr[i] =
            pow(-1.0, fnu_m + tmp) * (cyr_ptr[tmp] + I_IMPL_ * cyi_ptr[tmp]);
      }

      /* Free auxiliary pointers */
      free(--cyr_ptr);
      free(--cyi_ptr);
      free(--cwrkr_ptr);
      free(--cwrki_ptr);

    } else {

      /* Non-integer negative orders */

      /* Dynamic mem alloc of auxiliary arrays */
      double *cyr_ptr = (double *)malloc((n + 1) * sizeof(double));
      double *cyi_ptr = (double *)malloc((n + 1) * sizeof(double));
      double *cwrkr_ptr = (double *)malloc((n + 1) * sizeof(double));
      double *cwrki_ptr = (double *)malloc((n + 1) * sizeof(double));
      ++cyr_ptr;
      ++cyi_ptr;
      ++cwrkr_ptr;
      ++cwrki_ptr;

      /* Compute zbesy */
      double fnu_m = fabs(nu_m);
      slatec_zbesy_impl_(&x, &y, &fnu_m, &kode, &n, &cyr_ptr[0], &cyi_ptr[0],
                         &nz, &cwrkr_ptr[0], &cwrki_ptr[0], &ierr);
      slatec_flags_zbesy_impl_(ierr, nz);

      /* Dynamic mem alloc of auxiliary arrays */
      double *cjr_ptr = (double *)malloc((n + 1) * sizeof(double));
      double *cji_ptr = (double *)malloc((n + 1) * sizeof(double));
      ++cjr_ptr;
      ++cji_ptr;

      /* Compute zbesj */
      slatec_zbesj_impl_(&x, &y, &fnu_m, &kode, &n, &cjr_ptr[0], &cji_ptr[0],
                         &nz, &ierr);
      slatec_flags_zbesj_impl_(ierr, nz);

      /* Store in the array */
      for (int i = 0; i < n; i++) {
        int tmp1 = n - 1 - i;
        double tmp2 = (fnu_m + (double)tmp1) * M_PI;
        /* Eq. (5.5.4) of Ref. [2] */
        cyl_y_arr[i] = (cjr_ptr[tmp1] + I_IMPL_ * cji_ptr[tmp1]) * sin(tmp2) +
                       (cyr_ptr[tmp1] + I_IMPL_ * cyi_ptr[tmp1]) * cos(tmp2);
      }

      /* Free auxiliary pointers */
      free(--cyr_ptr);
      free(--cyi_ptr);
      free(--cwrkr_ptr);
      free(--cwrki_ptr);
      free(--cjr_ptr);
      free(--cji_ptr);
    }

  } else {

    /* Mixed negative and positive orders */

    int n_m = (int)floor(fabs(nu)) + 1;
    int n_p = n - n_m;
    double fnu_m = fabs(nu + (double)(n_m - 1));

    /* The negative orders */
    if (fabs(floor(fnu) - fnu) < DBL_EPSILON) {

      /* Negative integer orders, or order 0 */

      /* Dynamic mem alloc of auxiliary arrays */
      double *cyr_m_ptr = (double *)malloc((n_m + 1) * sizeof(double));
      double *cyi_m_ptr = (double *)malloc((n_m + 1) * sizeof(double));
      double *cwrkr_ptr = (double *)malloc((n_m + 1) * sizeof(double));
      double *cwrki_ptr = (double *)malloc((n_m + 1) * sizeof(double));
      ++cyr_m_ptr;
      ++cyi_m_ptr;
      ++cwrkr_ptr;
      ++cwrki_ptr;

      /* Compute zbesy */
      slatec_zbesy_impl_(&x, &y, &fnu_m, &kode, &n_m, &cyr_m_ptr[0],
                         &cyi_m_ptr[0], &nz, &cwrkr_ptr[0], &cwrki_ptr[0],
                         &ierr);
      slatec_flags_zbesy_impl_(ierr, nz);

      /* Store in the array */
      for (int i = 0; i < n_m; i++) {
        int tmp = n_m - 1 - i;
        cyl_y_arr[i] = pow(-1.0, fnu_m + tmp) *
                       (cyr_m_ptr[tmp] + I_IMPL_ * cyi_m_ptr[tmp]);
      }

      /* Free auxiliary pointers */
      free(--cyr_m_ptr);
      free(--cyi_m_ptr);
      free(--cwrkr_ptr);
      free(--cwrki_ptr);

    } else {

      /* Negative non-integer orders */

      /* Dynamic mem alloc of auxiliary arrays */
      double *cyr_m_ptr = (double *)malloc((n_m + 1) * sizeof(double));
      double *cyi_m_ptr = (double *)malloc((n_m + 1) * sizeof(double));
      double *cwrkr_ptr = (double *)malloc((n_m + 1) * sizeof(double));
      double *cwrki_ptr = (double *)malloc((n_m + 1) * sizeof(double));
      ++cyr_m_ptr;
      ++cyi_m_ptr;
      ++cwrkr_ptr;
      ++cwrki_ptr;

      /* Compute zbesy */
      slatec_zbesy_impl_(&x, &y, &fnu_m, &kode, &n_m, &cyr_m_ptr[0],
                         &cyi_m_ptr[0], &nz, &cwrkr_ptr[0], &cwrki_ptr[0],
                         &ierr);
      slatec_flags_zbesy_impl_(ierr, nz);

      /* Dynamic mem alloc of auxiliary arrays */
      double *cjr_m_ptr = (double *)malloc((n_m + 1) * sizeof(double));
      double *cji_m_ptr = (double *)malloc((n_m + 1) * sizeof(double));
      ++cjr_m_ptr;
      ++cji_m_ptr;

      /* Compute zbesj */
      slatec_zbesj_impl_(&x, &y, &fnu_m, &kode, &n_m, &cjr_m_ptr[0],
                         &cji_m_ptr[0], &nz, &ierr);
      slatec_flags_zbesj_impl_(ierr, nz);

      /* Store in the array */
      for (int i = 0; i < n_m; i++) {
        int tmp1 = n_m - 1 - i;
        double tmp2 = (fnu_m + (double)tmp1) * M_PI;
        /* Eq. (5.5.4) of Ref. [2] */
        cyl_y_arr[i] =
            (cjr_m_ptr[tmp1] + I_IMPL_ * cji_m_ptr[tmp1]) * sin(tmp2) +
            (cyr_m_ptr[tmp1] + I_IMPL_ * cyi_m_ptr[tmp1]) * cos(tmp2);
      }

      /* Free auxiliary pointers */
      free(--cyr_m_ptr);
      free(--cyi_m_ptr);
      free(--cwrkr_ptr);
      free(--cwrki_ptr);
      free(--cjr_m_ptr);
      free(--cji_m_ptr);
    }

    /* Positive orders */

    /* Auxiliary pointers */
    double *cyr_p_ptr = (double *)malloc((n_p + 1) * sizeof(double));
    double *cyi_p_ptr = (double *)malloc((n_p + 1) * sizeof(double));
    double *cwrkr_ptr = (double *)malloc((n_p + 1) * sizeof(double));
    double *cwrki_ptr = (double *)malloc((n_p + 1) * sizeof(double));
    ++cyr_p_ptr;
    ++cyi_p_ptr;
    ++cwrkr_ptr;
    ++cwrki_ptr;

    /* Compute zbesy */
    double fnu_p = nu + (double)n_m;
    slatec_zbesy_impl_(&x, &y, &fnu_p, &kode, &n_p, &cyr_p_ptr[0],
                       &cyi_p_ptr[0], &nz, &cwrkr_ptr[0], &cwrki_ptr[0], &ierr);
    slatec_flags_zbesy_impl_(ierr, nz);

    for (int i = n_m; i < n; i++) {
      int tmp1 = i - n_m;
      cyl_y_arr[i] = cyr_p_ptr[tmp1] + I_IMPL_ * cyi_p_ptr[tmp1];
    }

    /* Free auxiliary pointers */
    free(--cyr_p_ptr);
    free(--cyi_p_ptr);
    free(--cwrkr_ptr);
    free(--cwrki_ptr);
  }
}

#endif /* BESSEL_LIBRARY_CYL_Y_FULL_SEQ_IMPL_H */