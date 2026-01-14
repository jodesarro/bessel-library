/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/f2c_d_sign_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Function d_sign of the f2c library.
*/

#ifndef BESSEL_LIBRARY_F2C_D_SIGN_IMPL_H
#define BESSEL_LIBRARY_F2C_D_SIGN_IMPL_H

#include <math.h>

static inline double f2c_d_sign_impl_(double *x, double *y) {
    return (*y >= 0. ? fabs(*x) : -fabs(*x));
}

#endif /* BESSEL_LIBRARY_F2C_D_SIGN_IMPL_H */