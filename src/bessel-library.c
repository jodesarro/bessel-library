/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: src/bessel-library.c
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99 with C++ guards
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Wrapper for compiling the include/bessel-library.h file.
*/

/* Overwrite 'static inline' */
#if defined(_WIN32) || defined(_WIN64)
    #ifdef __cplusplus
        #define BESSEL_LIBRARY_STATIC_INLINE_IMPL_ __declspec(dllexport) extern "C"
    #else
        #define BESSEL_LIBRARY_STATIC_INLINE_IMPL_ __declspec(dllexport)
    #endif
#else
    #ifdef __cplusplus
        #define BESSEL_LIBRARY_STATIC_INLINE_IMPL_ extern "C"
    #else
        #define BESSEL_LIBRARY_STATIC_INLINE_IMPL_
    #endif
#endif

#include "../include/bessel-library.h"