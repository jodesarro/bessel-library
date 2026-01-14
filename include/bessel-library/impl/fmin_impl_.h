/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/fmin_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Returns the minimum between two real numbers of double type.
*/

#ifndef BESSEL_LIBRARY_FMIN_IMPL_H
#define BESSEL_LIBRARY_FMIN_IMPL_H

static inline double fmin_impl_(double x, double y) {
    return(x < y ? x : y);
}

#endif /* BESSEL_LIBRARY_FMIN_IMPL_H */