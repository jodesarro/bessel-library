/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/f2c_pow_dd_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt

    Description:
        Function pow_dd of the f2c library.
*/

#ifndef BESSEL_LIBRARY_F2C_POW_DD_IMPL_H
#define BESSEL_LIBRARY_F2C_POW_DD_IMPL_H

#include <math.h>

static inline double f2c_pow_dd_impl_(double *x, double *y) {
    return pow(*x, *y);
}

#endif /* BESSEL_LIBRARY_F2C_POW_DD_IMPL_H */