/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: src/bessel-library.c
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt

    Description:
        Wrapper for compiling the include/bessel-library.h
*/

#include <complex.h> /* for complex numbers */
typedef double complex tpdcomplex_impl_;
#include "../include/bessel-library/impl/airy_ai_impl_.h"
#include "../include/bessel-library/impl/airy_bi_impl_.h"
#include "../include/bessel-library/impl/cyl_h1_full_seq_impl_.h"
#include "../include/bessel-library/impl/cyl_h2_full_seq_impl_.h"
#include "../include/bessel-library/impl/cyl_i_full_seq_impl_.h"
#include "../include/bessel-library/impl/cyl_j_full_seq_impl_.h"
#include "../include/bessel-library/impl/cyl_k_full_seq_impl_.h"
#include "../include/bessel-library/impl/cyl_y_full_seq_impl_.h"

/*
    Returns the Airy function of the first kind and complex argument z, i.e.,
    Ai(z).

    Parameter:
    - z, complex argument of Ai(z).
*/
tpdcomplex_impl_ airy_ai(tpdcomplex_impl_ z) {
    return airy_ai_impl_(z, 0, 0);
}

/*
    Returns the first derivative of the Airy function of the first kind and
    complex argument z, i.e., dAi(z)/dz.

    Parameter:
    - z, complex argument of dAi(z)/dz.
*/
tpdcomplex_impl_ airy_ai_deriv(tpdcomplex_impl_ z) {
    return airy_ai_impl_(z, 1, 0);
}

/*
    Returns the scaled version of the Airy function of the first kind and
    complex argument z, i.e., Ai(z)*exp((2/3)*pow(z,3/2)).

    Parameter:
    - z, complex argument of Ai(z)*exp((2/3)*pow(z,3/2)).
*/
tpdcomplex_impl_ airy_ai_scal(tpdcomplex_impl_ z) {
    return airy_ai_impl_(z, 0, 1);
}

/*
    Returns the scaled version of the first derivative of the Airy function of
    the first kind and complex argument z, i.e.,
    (dAi(z)/dz)*exp((2/3)*pow(z,3/2)).

    Parameter:
    - z, complex argument of (dAi(z)/dz)*exp((2/3)*pow(z,3/2)).
*/
tpdcomplex_impl_ airy_ai_deriv_scal(tpdcomplex_impl_ z) {
    return airy_ai_impl_(z, 1, 1);
}

/*
    Returns the Airy function of the second kind and complex argument z, i.e.,
    Bi(z).

    Parameter:
    - z, complex argument of Bi(z).
*/
tpdcomplex_impl_ airy_bi(tpdcomplex_impl_ z) {
    return airy_bi_impl_(z, 0, 0);
}

/*
    Returns the first derivative of the Airy function of the second kind and
    complex argument z, i.e., dBi(z)/dz.

    Parameter:
    - z, complex argument of dBi(z)/dz.
*/
tpdcomplex_impl_ airy_bi_deriv(tpdcomplex_impl_ z) {
    return airy_bi_impl_(z, 1, 0);
}

/*
    Returns the scaled version of the Airy function of the second kind and
    complex argument z, i.e., Bi(z)*exp(-abs(real((2/3)*pow(z,3/2)))).

    Parameter:
    - z, complex argument of Bi(z)*exp(-abs(real((2/3)*pow(z,3/2)))).
*/
tpdcomplex_impl_ airy_bi_scal(tpdcomplex_impl_ z) {
    return airy_bi_impl_(z, 0, 1);
}

/*
    Returns the scaled version of the first derivative of the Airy function of
    the second kind and complex argument z, i.e.,
    (dBi(z)/dz)*exp(-abs(real((2/3)*pow(z,3/2)))).

    Parameter:
    - z, complex argument of (dBi(z)/dz)*exp(-abs(real((2/3)*pow(z,3/2)))).
*/
tpdcomplex_impl_ airy_bi_deriv_scal(tpdcomplex_impl_ z) {
    return airy_bi_impl_(z, 1, 1);
}

/*
    Returns the cylindrical Hankel function of the first kind, real order nu,
    and complex argument z, i.e., H1_nu(z).

    Parameters:
    - nu, real order of H1_nu(z).
    - z, complex argument of H1_nu(z).
*/
tpdcomplex_impl_ cyl_h1(double nu, tpdcomplex_impl_ z) {

    /* Array of one size */
    tpdcomplex_impl_ ch1[1];
    
    /* Compute cyl_h1_full_seq_impl_ */
    cyl_h1_full_seq_impl_(nu, 1, z, ch1, 0);

    /* Return */
    return ch1[0];
}

/*
    Returns the scaled version of the cylindrical Hankel function of the
    first kind, real order nu, and complex argument z, i.e.,
    H1_nu(z)*exp(-i*z).

    Parameters:
    - nu, real order of H1_nu(z)*exp(-i*z).
    - z, complex argument of H1_nu(z)*exp(-i*z).
*/
tpdcomplex_impl_ cyl_h1_scal(double nu,
    tpdcomplex_impl_ z) {

    /* Array of one size */
    tpdcomplex_impl_ ch1[1];
    
    /* Compute cyl_h1_full_seq_impl_ */
    cyl_h1_full_seq_impl_(nu, 1, z, ch1, 1);

    /* Return */
    return ch1[0];
}

/*
    Computes a n-sequency array of cylindrical Hankel functions of the first
    kind, real order nu, and complex argument z, i.e., {H1_nu(z),
    H1_(nu+1)(z), ..., H1_(nu+n-1)(z)}.

    Parameters:
    - nu, real order of H1_nu(z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_h1_arr array.
    - z, complex argument of H1_nu(z).
    - cyl_h1_arr, array of size n to output H1_nu(z) for the orders nu, nu+1,
    ..., nu+n-1
*/
void cyl_h1_seq(double nu, unsigned int n,
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_h1_arr) {

    cyl_h1_full_seq_impl_(nu, n, z, cyl_h1_arr, 0);
}

/*
    Computes a n-sequency array of the scaled version of the cylindrical
    Hankel functions of the first kind, real order nu, and complex argument z,
    i.e., {H1_nu(z)*exp(-i*z), H1_(nu+1)(z)*exp(-i*z), ...,
    H1_(nu+n-1)(z)*exp(-i*z)}.

    Parameters:
    - nu, real order of H1_nu(z)*exp(-i*z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_h2_scal_arr array.
    - z, complex argument of H1_nu(z)*exp(-i*z).
    - cyl_h2_scal_arr, array of size n to output  H1_nu(z)*exp(-i*z) for the
    orders nu, nu+1, ..., nu+n-1
*/
void cyl_h1_scal_seq(double nu, unsigned int n,
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_h2_scal_arr) {
    
    cyl_h1_full_seq_impl_(nu, n, z, cyl_h2_scal_arr, 1);
}

/*
    Returns the cylindrical Hankel function of the second kind, real order nu,
    and complex argument z, i.e., H2_nu(z).

    Parameters:
    - nu, real order of H2_nu(z).
    - z, complex argument of H2_nu(z).
*/
tpdcomplex_impl_ cyl_h2(double nu, tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ ch2[1];
    
    /* Compute cyl_h2_full_seq_impl_ */
    cyl_h2_full_seq_impl_(nu, 1, z, ch2, 0);

    /* Return */
    return ch2[0];
}

/*
    Returns the scaled version of the cylindrical Hankel function of the
    second kind, real order nu, and complex argument z, i.e.,
    H2_nu(z)*exp(i*z).

    Parameters:
    - nu, real order of H2_nu(z)*exp(i*z).
    - z, complex argument of H2_nu(z)*exp(i*z).
*/
tpdcomplex_impl_ cyl_h2_scal(double nu,
    tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ ch2[1];
    
    /* Compute cyl_h2_full_seq_impl_ */
    cyl_h2_full_seq_impl_(nu, 1, z, ch2, 1);

    /* Return */
    return ch2[0];
}

/*
    Computes a n-sequency array of cylindrical Hankel functions of the second
    kind, real order nu, and complex argument z, i.e., {H2_nu(z),
    H2_(nu+1)(z), ..., H2_(nu+n-1)(z)}.

    Parameters:
    - nu, real order of H2_nu(z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_h2_arr array.
    - z, complex argument of H2_nu(z).
    - cyl_h2_arr, array of size n to output H2_nu(z) for the orders nu, nu+1,
    ..., nu+n-1
*/
void cyl_h2_seq(double nu, unsigned int n,
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_h2_arr) {

    cyl_h2_full_seq_impl_(nu, n, z, cyl_h2_arr, 0);
}

/*
    Computes a n-sequency array of the scaled version of the cylindrical
    Hankel functions of the second kind, real order nu, and complex argument
    z, i.e., {H1_nu(z)*exp(i*z), H1_(nu+1)(z)*exp(i*z), ...,
    H1_(nu+n-1)(z)*exp(i*z)}.

    Parameters:
    - nu, real order of H2_nu(z)*exp(i*z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_h2_scal_arr array.
    - z, complex argument of H2_nu(z)*exp(i*z).
    - cyl_h2_scal_arr, array of size n to output H2_nu(z)*exp(i*z) for the
    orders nu, nu+1, ..., nu+n-1
*/
void cyl_h2_scal_seq(double nu, unsigned int n,
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_h2_scal_arr) {
    
    cyl_h2_full_seq_impl_(nu, n, z, cyl_h2_scal_arr, 1);
}

/*
    Returns the modified cylindrical Bessel function of the first kind, real
    order nu, and complex argument z, i.e., I_nu(z).

    Parameters:
    - nu, real order of I_nu(z).
    - z, complex argument of I_nu(z).
*/
tpdcomplex_impl_ cyl_i(double nu, tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ ci[1];
    
    /* Compute cyl_i_full_seq_impl_ */
    cyl_i_full_seq_impl_(nu, 1, z, ci, 0);

    /* Return */
    return ci[0];
}

/*
    Returns the scaled version of the modified cylindrical Bessel function of
    the first kind, real order nu, and complex argument z, i.e.,
    I_nu(z)*exp(-abs(real(z))).

    Parameters:
    - nu, real order of I_nu(z)*exp(-abs(real(z))).
    - z, complex argument of I_nu(z)*exp(-abs(real(z))).
*/
tpdcomplex_impl_ cyl_i_scal(double nu,
    tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ ci[1];
    
    /* Compute cyl_i_full_seq_impl_ */
    cyl_i_full_seq_impl_(nu, 1, z, ci, 1);

    /* Return */
    return ci[0];
}

/*
    Computes a n-sequency array of modified cylindrical Bessel functions of
    the first kind, real order nu, and complex argument z, i.e.,
    {I_nu(z), I_(nu+1)(z), ..., I_(nu+n-1)(z)}.

    Parameters:
    - nu, real order of I_nu(z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_i_arr array.
    - z, complex argument of I_nu(z).
    - cyl_i_arr, array of size n to output I_nu(z) for the orders nu, nu+1,
    ..., nu+n-1
*/
static inline
void cyl_i_seq(double nu, unsigned int n, tpdcomplex_impl_ z,
    tpdcomplex_impl_ *cyl_i_arr) {

    cyl_i_full_seq_impl_(nu, n, z, cyl_i_arr, 0);
}

/*
    Computes a n-sequency array of the scaled version of the modified
    cylindrical Bessel functions of the first kind, real order nu, and complex
    argument z, i.e., {I_nu(z)*exp(-abs(real(z))),
    I_(nu+1)(z)*exp(-abs(real(z))), ..., I_(nu+n-1)(z)*exp(-abs(real(z)))}.

    Parameters:
    - nu, real order of I_nu(z)*exp(-abs(real(z))).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_i_scal_arr array.
    - z, complex argument of I_nu(z)*exp(-abs(real(z))).
    - cyl_i_scal_arr, array of size n to output I_nu(z)*exp(-abs(real(z))) for
    the orders nu, nu+1, ..., nu+n-1
*/
void cyl_i_scal_seq(double nu, unsigned int n,
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_i_scal_arr) {
    
    cyl_i_full_seq_impl_(nu, n, z, cyl_i_scal_arr, 1);
}

/*
    Returns the modified cylindrical Bessel function of the second kind, real
    order nu, and complex argument z, i.e., K_nu(z).

    Parameters:
    - nu, real order of K_nu(z).
    - z, complex argument of K_nu(z).
*/
tpdcomplex_impl_ cyl_k(double nu, tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ ck[1];
    
    /* Compute cyl_k_full_seq_impl_ */
    cyl_k_full_seq_impl_(nu, 1, z, ck, 0);

    /* Return */
    return ck[0];
}

/*
    Returns the scaled version of the modified cylindrical Bessel function of
    the second kind, real order nu, and complex argument z, i.e.,
    K_nu(z)*exp(z).

    Parameters:
    - nu, real order of K_nu(z)*exp(z).
    - z, complex argument of K_nu(z)*exp(z).
*/
tpdcomplex_impl_ cyl_k_scal(double nu,
    tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ ck[1];
    
    /* Compute cyl_k_full_seq_impl_ */
    cyl_k_full_seq_impl_(nu, 1, z, ck, 1);

    /* Return */
    return ck[0];
}

/*
    Computes a n-sequency array of modified cylindrical Bessel functions of
    the second kind, real order nu, and complex argument z, i.e.,
    {K_nu(z), Y_(nu+1)(z), ..., Y_(nu+n-1)(z)}.

    Parameters:
    - nu, real order of K_nu(z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_k_arr array.
    - z, complex argument of K_nu(z).
    - cyl_k_arr, array of size n to output K_nu(z) for the orders nu,
    nu+1, ..., nu+n-1
*/
void cyl_k_seq(double nu, unsigned int n, tpdcomplex_impl_ z,
    tpdcomplex_impl_ *cyl_k_arr) {

    cyl_k_full_seq_impl_(nu, n, z, cyl_k_arr, 0);
}

/*
    Computes a n-sequency array of the scaled version of the modified
    cylindrical Bessel functions of the second kind, real order nu, and
    complex argument z, i.e., {K_nu(z)*exp(z), Y_(nu+1)(z)*exp(z),
    ..., Y_(nu+n-1)(z)*exp(z)}.

    Parameters:
    - nu, real order of K_nu(z)*exp(z).
    - n, number of elements in the sequence for computing the orders nu,
    nu+1, ..., nu+n-1 It is also the size of the cyl_k_scaled_arr array.
    - z, complex argument of K_nu(z)*exp(z).
    - cyl_k_scaled_arr, array of size n to output K_nu(z)*exp(z) for the
    orders nu, nu+1, ..., nu+n-1
*/
void cyl_k_scal_seq(double nu, unsigned int n,
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_k_scaled_arr) {
    
    cyl_k_full_seq_impl_(nu, n, z, cyl_k_scaled_arr, 1);
}

/*
    Returns the cylindrical Bessel function of the second kind, real order nu,
    and complex argument z, i.e., Y_nu(z).

    Parameters:
    - nu, real order of Y_nu(z).
    - z, complex argument of Y_nu(z).
*/
tpdcomplex_impl_ cyl_y(double nu, tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ cy[1];
    
    /* Compute cyl_y_full_seq_impl_ */
    cyl_y_full_seq_impl_(nu, 1, z, cy, 0);

    /* Return */
    return cy[0];
}

/*
    Returns the scaled version of the cylindrical Bessel function of the
    second kind, real order nu, and complex argument z, i.e.,
    Y_nu(z)*exp(-abs(imag(z))).

    Parameters:
    - nu, real order of Y_nu(z)*exp(-abs(imag(z))).
    - z, complex argument of Y_nu(z)*exp(-abs(imag(z))).
*/
tpdcomplex_impl_ cyl_y_scal(double nu,
    tpdcomplex_impl_ z) {
    
    /* Array of one size */
    tpdcomplex_impl_ cy[1];
    
    /* Compute cyl_y_full_seq_impl_ */
    cyl_y_full_seq_impl_(nu, 1, z, cy, 1);

    /* Return */
    return cy[0];
}

/*
    Computes a n-sequency array of cylindrical Bessel functions of the second
    kind, real order nu, and complex argument z, i.e.,
    {Y_nu(z), Y_(nu+1)(z), ..., Y_(nu+n-1)(z)}.

    Parameters:
    - nu, real order of Y_nu(z).
    - n, number of elements in the sequence for computing the orders nu, nu+1,
    ..., nu+n-1 It is also the size of the cyl_y_arr array.
    - z, complex argument of Y_nu(z).
    - cyl_y_arr, array of size n to output Y_nu(z) for the orders nu,
    nu+1, ..., nu+n-1
*/
void cyl_y_seq(double nu, unsigned int n, tpdcomplex_impl_ z,
    tpdcomplex_impl_ *cyl_y_arr) {

    cyl_y_full_seq_impl_(nu, n, z, cyl_y_arr, 0);
}

/*
    Computes a n-sequency array of the scaled version of the cylindrical
    Bessel functions of the second kind, real order nu, and complex argument
    z, i.e., {Y_nu(z)*exp(-abs(imag(z))), Y_(nu+1)(z)*exp(-abs(imag(z))),
    ..., Y_(nu+n-1)(z)*exp(-abs(imag(z)))}.

    Parameters:
    - nu, real order of Y_nu(z)*exp(-abs(imag(z))).
    - n, number of elements in the sequence for computing the orders nu,
    nu+1, ..., nu+n-1 It is also the size of the cyl_y_scaled_arr array.
    - z, complex argument of Y_nu(z)*exp(-abs(imag(z))).
    - cyl_y_scaled_arr, array of size n to output Y_nu(z)*exp(-abs(imag(z)))
    for the orders nu, nu+1, ..., nu+n-1
*/
void cyl_y_scal_seq(double nu, unsigned int n,
    tpdcomplex_impl_ z, tpdcomplex_impl_ *cyl_y_scaled_arr) {
    
    cyl_y_full_seq_impl_(nu, n, z, cyl_y_scaled_arr, 1);
}