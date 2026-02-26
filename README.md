# Bessel Library: A C library with routines for computing Bessel functions

This project provides a C library
([compatible with C++](#compatibility-with-c)) designed for the
computation of Bessel functions and related special functions for real orders
and complex arguments.

Detailed discussions on Bessel functions and their numerical evaluation may be
found in the works listed in the [References](#references) section.

## Available features

### $J_\nu(z)$—Cylindrical Bessel function of the first kind

<details>
  <summary>
    <code><b>cyl_j(nu, z)</b></code>
  </summary>

  - **Description:** Returns the cylindrical Bessel function of the first
  kind, real order $\nu$, and complex argument $z$, i.e., $J_\nu(z)$.
  - **Parameters:**
    - `nu`, real order of $J_\nu(z)$.
    - `z`, complex argument of $J_\nu(z)$.
  - **Implementation:** Similar to the `cyl_j_seq()` function.
</details>

<details>
  <summary>
    <code><b>cyl_j_scal(nu, z)</b></code> 
  </summary>

  - **Description:** Returns a scaled version of the cylindrical Bessel
  function of the first kind, real order $\nu$, and complex argument $z$,
  i.e., $J_\nu(z) \ e^{-|\mathrm{Im}(z)|}$.
  - **Parameters:**
    - `nu`, real order of $J_\nu(z) \ e^{-|\mathrm{Im}(z)|}$.
    - `z`, complex argument of $J_\nu(z) \ e^{-|\mathrm{Im}(z)|}$.
  - **Implementation:** Similar to the `cyl_j_seq()` function.
</details>

<details>
  <summary>
    <code><b>cyl_j_seq(nu, n, z, cyl_j_arr)</b></code>
  </summary>

  - **Description:** Computes a $n$-sequency array of cylindrical Bessel
  functions of the first kind, real order $\nu$, and complex argument $z$,
  i.e., { $J_\nu(z)$, $J_{\nu+1}(z)$, ..., $J_{\nu+n-1}(z)$ }.
  - **Parameters:**
    - `nu`, real order of $J_\nu(z)$.
    - `n`, number $n$ of elements in the sequence for computing the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$. It is also the size of the `cyl_j_arr`
    array.
    - `z`, complex argument of $J_\nu(z)$.
    - `cyl_j_arr`, array of size $n$ to output $J_\nu(z)$ for the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$.
  - **Implementation:** In general, the implementation is based on the D. E.
  Amos Fortran 77 routines from the Slatec library [[3]](#ref3). Such
  Fortran routines, and all their dependencies, were carefully translated to
  C. Negative orders are handled by Eqs. (5.4.2) and (5.5.4) of
  Ref. [[2]](#ref2) for, respectively, $\nu \in \mathtt{Z}$ and
  $\nu \notin \mathtt{Z}$; in the latter case, it yields $\infty+i\infty$ when
  $|z|=0$.
</details>

<details>
  <summary>
    <code><b>cyl_j_scal_seq(nu, n, z, cyl_j_scal_arr)</b></code> 
  </summary>

  - **Description:** Computes a $n$-sequency array of scaled versions of
  cylindrical Bessel functions of the first kind, real order $\nu$, and
  complex argument $z$, i.e., { $J_\nu(z) \ e^{-|\mathrm{Im}(z)|}$,
  $J_{\nu+1}(z) \ e^{-|\mathrm{Im}(z)|}$, ...,
  $J_{\nu+n-1}(z) \ e^{-|\mathrm{Im}(z)|}$ }.
  - **Parameters:**
    - `nu`, real order of $J_\nu(z) \ e^{-|\mathrm{Im}(z)|}$.
    - `n`, number $n$ of elements in the sequence for computing the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$. It is also the size of the
    `cyl_j_scal_arr` array.
    - `z`, complex argument of $J_\nu(z) \ e^{-|\mathrm{Im}(z)|}$.
    - `cyl_j_scal_arr`, array of size $n$ to output
    $J_\nu(z) \ e^{-|\mathrm{Im}(z)|}$ for the orders $\nu$, $\nu+1$, ...,
    $\nu+n-1$.
  - **Implementation:** Similar to the `cyl_j_seq()` function.
</details>

### $Y_\nu(z)$—Cylindrical Bessel function of the second kind

<details>
  <summary>
    <b><code>cyl_y(nu, z)</code></b>
  </summary>

  - **Description:** Returns the cylindrical Bessel function of the second
  kind, real order $\nu$, and complex argument $z$, i.e., $Y_\nu(z)$.
  - **Parameters:**
    - `nu`, real order of $Y_\nu(z)$.
    - `z`, complex argument of $Y_\nu(z)$.
  - **Implementation:** Similar to the `cyl_y_seq()` function.
</details>

<details>
  <summary>
    <b><code>cyl_y_scal(nu, z)</code></b>
  </summary>

  - **Description:** Returns a scaled version of the cylindrical Bessel
  function of the second kind, real order $\nu$, and complex argument $z$,
  i.e., $Y_\nu(z) \ e^{-|\mathrm{Im}(z)|}$.
  - **Parameters:**
    - `nu`, real order of $Y_\nu(z) \ e^{-|\mathrm{Im}(z)|}$.
    - `z`, complex argument of $Y_\nu(z) \ e^{-|\mathrm{Im}(z)|}$.
  - **Implementation:** Similar to the `cyl_y_seq()` function.
</details>

<details>
  <summary>
    <b><code>cyl_y_seq(nu, n, z, cyl_y_arr)</code></b>
  </summary>

  - **Description:** Computes a $n$-sequency array of cylindrical Bessel
  functions of the second kind, real order $\nu$, and complex argument $z$,
  i.e., { $Y_\nu(z)$, $Y_{\nu+1}(z)$, ..., $Y_{\nu+n-1}(z)$ }.
  - **Parameters:**
    - `nu`, real order of $Y_\nu(z)$.
    - `n`, number $n$ of elements in the sequence for computing the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$. It is also the size of the `cyl_y_arr`
    array.
    - `z`, complex argument of $Y_\nu(z)$.
    - `cyl_y_arr`, array of size $n$ to output $Y_\nu(z)$ for the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$.
  - **Implementation:** In general, the implementation is based on the D. E.
  Amos Fortran 77 routines from the Slatec library [[3]](#ref3). Such
  Fortran routines, and all their dependencies, were carefully translated to
  C. Negative orders are handled by Eqs. (5.4.2) and (5.5.4) of
  Ref. [[2]](#ref2) for, respectively, $\nu \in \mathtt{Z}$ and
  $\nu \notin \mathtt{Z}$. When $|z|=0$, it yields $-\infty$ if $\nu=0$, or
  $\infty+i\infty$ otherwise.
</details>

<details>
  <summary>
    <b><code>cyl_y_scal_seq(nu, n, z, cyl_y_scal_arr)</code></b>
  </summary>

  - **Description:** Computes a $n$-sequency array of scaled versions of
  cylindrical Bessel functions of the second kind, real order $\nu$, and
  complex argument $z$, i.e., { $Y_\nu(z) \ e^{-|\mathrm{Im}(z)|}$,
  $Y_{\nu+1}(z) \ e^{-|\mathrm{Im}(z)|}$, ...,
  $Y_{\nu+n-1}(z) \ e^{-|\mathrm{Im}(z)|}$ }.
  - **Parameters:**
    - `nu`, real order of $Y_\nu(z) \ e^{-|\mathrm{Im}(z)|}$.
    - `n`, number $n$ of elements in the sequence for computing the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$. It is also the size of the
    `cyl_y_scal_arr` array.
    - `z`, complex argument of $Y_\nu(z) \ e^{-|\mathrm{Im}(z)|}$.
    - `cyl_j_scal_arr`, array of size $n$ to output
    $Y_\nu(z) \ e^{-|\mathrm{Im}(z)|}$ for the orders $\nu$, $\nu+1$, ...,
    $\nu+n-1$.
  - **Implementation:** Similar to the `cyl_y_seq()` function.
</details>

### $H_\nu^{(1)}(z)$—Cylindrical Hankel function of the first kind

<details>
  <summary>
    <b><code>cyl_h1(nu, z)</code></b>
  </summary>

  - **Description:** Returns the cylindrical Hankel function of the first
  kind, real order $\nu$, and complex argument $z$, i.e., $H_\nu^{(1)}(z)$.
  - **Parameters:**
    - `nu`, real order of $H_\nu^{(1)}(z)$.
    - `z`, complex argument of $H_\nu^{(1)}(z)$.
  - **Implementation:** Similar to the `cyl_h1_seq()` function.
</details>

<details>
  <summary>
    <b><code>cyl_h1_scal(nu, z)</code></b>
  </summary>

  - **Description:** Returns a scaled version of the cylindrical Hankel
  function of the first kind, real order $\nu$, and complex argument $z$,
  i.e., $H_\nu^{(1)}(z) \ e^{-iz}$.
  - **Parameters:**
    - `nu`, real order of $H_\nu^{(1)}(z) \ e^{-iz}$.
    - `z`, complex argument of $H_\nu^{(1)}(z) \ e^{-iz}$.
  - **Implementation:** Similar to the `cyl_h1_seq()` function.
</details>

<details>
  <summary>
    <b><code>cyl_h1_seq(nu, n, z, cyl_h1_arr)</code></b>
  </summary>

  - **Description:** Computes a $n$-sequency array of cylindrical Hankel
  functions of the first kind, real order $\nu$, and complex argument $z$,
  i.e., { $H_\nu^{(1)}(z)$, $H_{\nu+1}^{(1)}(z)$, ...,
  $H_{\nu+n-1}^{(1)}(z)$ }.
  - **Parameters:**
    - `nu`, real order of $H_\nu^{(1)}(z)$.
    - `n`, number $n$ of elements in the sequence for computing the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$. It is also the size of the `cyl_h1_arr`
    array.
    - `z`, complex argument of $H_\nu^{(1)}(z)$.
    - `cyl_h1_arr`, array of size $n$ to output $H_\nu^{(1)}(z)$ for the
    orders $\nu$, $\nu+1$, ..., $\nu+n-1$.
  - **Implementation:** In general, the implementation is based on the D. E.
  Amos Fortran 77 routines from the Slatec library [[3]](#ref3). Such
  Fortran routines, and all their dependencies, were carefully translated to
  C. Negative orders are handled by Eq. (9.1.6) of Ref. [[1]](#ref1). It
  yields $\infty+i\infty$ when $|z|=0$.
</details>

<details>
  <summary>
    <b><code>cyl_h1_scal_seq(nu, n, z, cyl_h1_scal_arr)</code></b>
  </summary>

  - **Description:** Computes a $n$-sequency array of scaled versions of
  cylindrical Hankel functions of the first kind, real order $\nu$, and
  complex argument $z$, i.e., { $H_\nu^{(1)}(z) \ e^{-iz}$,
  $H_{\nu+1}^{(1)}(z) \ e^{-iz}$, ..., $H_{\nu+n-1}^{(1)}(z) \ e^{-iz}$ }.
  - **Parameters:**
    - `nu`, real order of $H_\nu^{(1)}(z) \ e^{-iz}$.
    - `n`, number $n$ of elements in the sequence for computing the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$. It is also the size of the
    `cyl_h1_scal_arr` array.
    - `z`, complex argument of $H_\nu^{(1)}(z) \ e^{-iz}$.
    - `cyl_h1_scal_arr`, array of size $n$ to output
    $H_\nu^{(1)}(z) \ e^{-iz}$
    for the orders $\nu$, $\nu+1$, ..., $\nu+n-1$.
  - **Implementation:** Similar to the `cyl_h1_seq()` function.
</details>

### $H_\nu^{(2)}(z)$—Cylindrical Hankel function of the second kind

<details>
  <summary>
    <b><code>cyl_h2(nu, z)</code></b>
  </summary>

  - **Description:** Returns the cylindrical Hankel function of the second
  kind, real order $\nu$, and complex argument $z$, i.e., $H_\nu^{(2)}(z)$.
  - **Parameters:**
    - `nu`, real order of $H_\nu^{(2)}(z)$.
    - `z`, complex argument of $H_\nu^{(2)}(z)$.
  - **Implementation:** Similar to the `cyl_h2_seq()` function.
</details>

<details>
  <summary>
    <b><code>cyl_h2_scal(nu, z)</code></b>
  </summary>

  - **Description:** Returns a scaled version of the cylindrical Hankel
  function of the second kind, real order $\nu$, and complex argument $z$,
  i.e., $H_\nu^{(2)}(z) \ e^{iz}$.
  - **Parameters:**
    - `nu`, real order of $H_\nu^{(2)}(z) \ e^{iz}$.
    - `z`, complex argument of $H_\nu^{(2)}(z) \ e^{iz}$.
  - **Implementation:** Similar to the `cyl_h2_seq()` function.
</details>

<details>
  <summary>
    <b><code>cyl_h2_seq(nu, n, z, cyl_h2_arr)</code></b>
  </summary>

  - **Description:** Computes a $n$-sequency array of cylindrical Hankel
  functions of the second kind, real order $\nu$, and complex argument $z$,
  i.e., { $H_\nu^{(2)}(z)$, $H_{\nu+1}^{(2)}(z)$, ...,
  $H_{\nu+n-1}^{(2)}(z)$ }.
  - **Parameters:**
    - `nu`, real order of $H_\nu^{(2)}(z)$.
    - `n`, number $n$ of elements in the sequence for computing the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$. It is also the size of the `cyl_h2_arr`
    array.
    - `z`, complex argument of $H_\nu^{(2)}(z)$.
    - `cyl_h2_arr`, array of size $n$ to output $H_\nu^{(2)}(z)$ for the
    orders $\nu$, $\nu+1$, ..., $\nu+n-1$.
  - **Implementation:** In general, the implementation is based on the D. E.
  Amos Fortran 77 routines from the Slatec library [[3]](#ref3). Such
  Fortran routines, and all their dependencies, were carefully translated to
  C. Negative orders are handled by Eq. (9.1.6) of Ref. [[1]](#ref1). It
  yields $\infty+i\infty$ when $|z|=0$.
</details>

<details>
  <summary>
    <b><code>cyl_h2_scal_seq(nu, n, z, cyl_h2_scal_arr)</code></b>
  </summary>

  - **Description:** Computes a $n$-sequency array of scaled versions of
  cylindrical Hankel functions of the second kind, real order $\nu$, and
  complex argument $z$, i.e., { $H_\nu^{(2)}(z) \ e^{iz}$,
  $H_{\nu+1}^{(2)}(z) \ e^{iz}$, ..., $H_{\nu+n-1}^{(2)}(z) \ e^{iz}$ }.
  - **Parameters:**
    - `nu`, real order of $H_\nu^{(2)}(z) \ e^{iz}$.
    - `n`, number $n$ of elements in the sequence for computing the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$. It is also the size of the
    `cyl_h2_scal_arr` array.
    - `z`, complex argument of $H_\nu^{(2)}(z) \ e^{iz}$.
    - `cyl_h2_scal_arr`, array of size $n$ to output $H_\nu^{(2)}(z) \ e^{iz}$
    for the orders $\nu$, $\nu+1$, ..., $\nu+n-1$.
  - **Implementation:** Similar to the `cyl_h2_seq()` function.
</details>

### $I_\nu(z)$—Modified cylindrical Bessel function of the first kind

<details>
  <summary>
    <b><code>cyl_i(nu, z)</code></b>
  </summary>

  - **Description:** Returns the modified cylindrical Bessel function of the
  first kind, real order $\nu$, and complex argument $z$, i.e., $I_\nu(z)$.
  - **Parameters:**
    - `nu`, real order of $I_\nu(z)$.
    - `z`, complex argument of $I_\nu(z)$.
  - **Implementation:** Similar to the `cyl_i_seq()` function.
</details>

<details>
  <summary>
    <b><code>cyl_i_scal(nu, z)</code></b>
  </summary>

  - **Description:** Returns a scaled version of the modified cylindrical
  Bessel function of the first kind, real order $\nu$, and complex argument
  $z$, i.e., $I_\nu(z) \ e^{-|\mathrm{Re}(z)|}$.
  - **Parameters:**
    - `nu`, real order of $I_\nu(z) \ e^{-|\mathrm{Re}(z)|}$.
    - `z`, complex argument of $I_\nu(z) \ e^{-|\mathrm{Re}(z)|}$.
  - **Implementation:** Similar to the `cyl_i_seq()` function.
</details>

<details>
  <summary>
    <b><code>cyl_i_seq(nu, n, z, cyl_i_arr)</code></b>
  </summary>

  - **Description:** Computes a $n$-sequency array of modified cylindrical
  Bessel functions of the first kind, real order $\nu$, and complex argument
  $z$, i.e., { $I_\nu(z)$, $I_{\nu+1}(z)$, ..., $I_{\nu+n-1}(z)$ }.
  - **Parameters:**
    - `nu`, real order of $I_\nu(z)$.
    - `n`, number $n$ of elements in the sequence for computing the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$. It is also the size of the `cyl_i_arr`
    array.
    - `z`, complex argument of $I_\nu(z)$.
    - `cyl_i_arr`, array of size $n$ to output $I_\nu(z)$ for the
    orders $\nu$, $\nu+1$, ..., $\nu+n-1$.
  - **Implementation:** In general, the implementation is based on the D. E.
  Amos Fortran 77 routines from the Slatec library [[3]](#ref3). Such
  Fortran routines, and all their dependencies, were carefully translated to
  C. Negative orders are handled by Eqs. (6.1.5) and (6.5.4) of
  Ref. [[2]](#ref2) for, respectively, $\nu \in \mathtt{Z}$ and
  $\nu \notin \mathtt{Z}$; in the latter case, it yields $\infty+i\infty$ when
  $|z|=0$.
</details>

<details>
  <summary>
    <b><code>cyl_i_scal_seq(nu, n, z, cyl_i_scal_arr)</code></b>
  </summary>

  - **Description:** Computes a $n$-sequency array of scaled versions of
  modified cylindrical Bessel functions of the first kind, real order
  $\nu$, and complex argument $z$, i.e., { $I_\nu(z) \ e^{-|\mathrm{Re}(z)|}$,
  $I_{\nu+1}(z) \ e^{-|\mathrm{Re}(z)|}$, ...,
  $I_{\nu+n-1}(z) \ e^{-|\mathrm{Re}(z)|}$ }.
  - **Parameters:**
    - `nu`, real order of $I_\nu(z) \ e^{-|\mathrm{Re}(z)|}$.
    - `n`, number $n$ of elements in the sequence for computing the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$. It is also the size of the
    `cyl_i_scal_arr` array.
    - `z`, complex argument of $I_\nu(z) \ e^{-|\mathrm{Re}(z)|}$.
    - `cyl_i_scal_arr`, array of size $n$ to output
    $I_\nu(z) \ e^{-|\mathrm{Re}(z)|}$ for the orders $\nu$, $\nu+1$, ...,
    $\nu+n-1$.
  - **Implementation:** Similar to the `cyl_i_seq()` function.
</details>

### $K_\nu(z)$—Modified cylindrical Bessel function of the second kind

<details>
  <summary>
    <b><code>cyl_k(nu, z)</code></b>
  </summary>

  - **Description:** Returns the modified cylindrical Bessel function of the
  second kind, real order $\nu$, and complex argument $z$, i.e., $K_\nu(z)$.
  - **Parameters:**
    - `nu`, real order of $K_\nu(z)$.
    - `z`, complex argument of $K_\nu(z)$.
  - **Implementation:** Similar to the `cyl_k_seq()` function.
</details>

<details>
  <summary>
    <b><code>cyl_k_scal(nu, z)</code></b>
  </summary>

  - **Description:** Returns a scaled version of the modified cylindrical
  Bessel function of the second kind, real order $\nu$, and complex argument
  $z$, i.e., $K_\nu(z) \ e^{z}$.
  - **Parameters:**
    - `nu`, real order of $K_\nu(z) \ e^{z}$.
    - `z`, complex argument of $K_\nu(z) \ e^{z}$.
  - **Implementation:** Similar to the `cyl_k_seq()` function.
</details>

<details>
  <summary>
    <b><code>cyl_k_seq(nu, n, z, cyl_k_arr)</code></b>
  </summary>

  - **Description:** Computes a $n$-sequency array of modified cylindrical
  Bessel functions of the second kind, real order $\nu$, and complex argument
  $z$, i.e., { $K_\nu(z)$, $K_{\nu+1}(z)$, ..., $K_{\nu+n-1}(z)$ }.
  - **Parameters:**
    - `nu`, real order of $K_\nu(z)$.
    - `n`, number $n$ of elements in the sequence for computing the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$. It is also the size of the `cyl_k_arr`
    array.
    - `z`, complex argument of $K_\nu(z)$.
    - `cyl_k_arr`, array of size $n$ to output $K_\nu(z)$ for the
    orders $\nu$, $\nu+1$, ..., $\nu+n-1$.
  - **Implementation:** In general, the implementation is based on the D. E.
  Amos Fortran 77 routines from the Slatec library [[3]](#ref3). Such
  Fortran routines, and all their dependencies, were carefully translated to
  C. Negative orders are handled by Eqs. (6.5.5) of Ref. [[2]](#ref2).
  When $|z|=0$, it yields $\infty$ if $\nu=0$, or $\infty+i\infty$ otherwise.
</details>

<details>
  <summary>
    <b><code>cyl_k_scal_seq(nu, n, z, cyl_k_scal_arr)</code></b>
  </summary>

  - **Description:** Computes a $n$-sequency array of scaled versions of
  modified cylindrical Bessel functions of the second kind, real order
  $\nu$, and complex argument $z$, i.e., { $K_\nu(z) \ e^{z}$,
  $K_{\nu+1}(z) \ e^{z}$, ..., $K_{\nu+n-1}(z) \ e^{z}$ }.
  - **Parameters:**
    - `nu`, real order of $K_\nu(z) \ e^{z}$.
    - `n`, number $n$ of elements in the sequence for computing the orders
    $\nu$, $\nu+1$, ..., $\nu+n-1$. It is also the size of the
    `cyl_k_scal_arr` array.
    - `z`, complex argument of $K_\nu(z) \ e^{z}$.
    - `cyl_k_scal_arr`, array of size $n$ to output
    $K_\nu(z) \ e^{z}$ for the orders $\nu$, $\nu+1$, ...,
    $\nu+n-1$.
  - **Implementation:** Similar to the `cyl_k_seq()` function.
</details>

### $Ai(z)$—Airy function of the first kind

<details>
  <summary>
    <b><code>airy_ai(z)</code></b>
  </summary>

  - **Description:** Returns the Airy function of the first kind and complex
  argument $z$, i.e., $Ai(z)$.
  - **Parameter:**
    - `z`, complex argument of $Ai(z)$.
  - **Implementation:** In general, the implementation is based on the D. E.
  Amos Fortran 77 routines from the Slatec library [[3]](#ref3). Such
  Fortran routines, and all their dependencies, were carefully translated to
  C.
</details>

<details>
  <summary>
    <b><code>airy_ai_diff(z)</code></b>
  </summary>

  - **Description:** Returns the first derivative of the Airy function of the
  first kind and complex argument $z$, i.e., $\mathrm{d}Ai(z)/\mathrm{d}z$.
  - **Parameter:**
    - `z`, complex argument of $\mathrm{d}Ai(z)/\mathrm{d}z$.
  - **Implementation:** Similar to the `airy_ai()` function.
</details>

<details>
  <summary>
    <b><code>airy_ai_scal(z)</code></b>
  </summary>

  - **Description:** Returns a scaled version of the Airy function of the
  first kind and complex argument $z$, i.e., $Ai(z) \ e^{(2/3)z^{3/2}}$.
  - **Parameter:**
    - `z`, complex argument of $Ai(z) \ e^{(2/3)z^{3/2}}$.
  - **Implementation:** Similar to the `airy_ai()` function.
</details>

<details>
  <summary>
    <b><code>airy_ai_diff_scal(z)</code></b>
  </summary>

  - **Description:** Returns a scaled version of the first derivative of the
  Airy function of the first kind and complex argument $z$, i.e.,
  $[\mathrm{d}Ai(z)/\mathrm{d}z] \ e^{(2/3)z^{3/2}}$.
  - **Parameter:**
    - `z`, complex argument of
    $[\mathrm{d}Ai(z)/\mathrm{d}z] \ e^{(2/3)z^{3/2}}$.
  - **Implementation:** Similar to the `airy_ai()` function.
</details>

### $Bi(z)$—Airy function of the second kind

<details>
  <summary>
    <b><code>airy_bi(z)</code></b>
  </summary>

  - **Description:** Returns the Airy function of the second kind and complex
  argument $z$, i.e., $Bi(z)$.
  - **Parameter:**
    - `z`, complex argument of $Bi(z)$.
  - **Implementation:** In general, the implementation is based on the D. E.
  Amos Fortran 77 routines from the Slatec library [[3]](#ref3). Such
  Fortran routines, and all their dependencies, were carefully translated to
  C.
</details>

<details>
  <summary>
    <b><code>airy_bi_diff(z)</code></b>
  </summary>

  - **Description:** Returns the first derivative of the Airy function of the
  second kind and complex argument $z$, i.e., $\mathrm{d}Bi(z)/\mathrm{d}z$.
  - **Parameter:**
    - `z`, complex argument of $\mathrm{d}Bi(z)/\mathrm{d}z$.
  - **Implementation:** Similar to the `airy_bi()` function.
</details>

<details>
  <summary>
    <b><code>airy_ai_scal(z)</code></b>
  </summary>

  - **Description:** Returns a scaled version of the Airy function of the
  second kind and complex argument $z$, i.e.,
  $Bi(z) \ e^{-|\mathrm{Re}[(2/3)z^{3/2}]|}$.
  - **Parameter:**
    - `z`, complex argument of $Bi(z) \ e^{-|\mathrm{Re}[(2/3)z^{3/2}]|}$.
  - **Implementation:** Similar to the `airy_bi()` function.
</details>

<details>
  <summary>
    <b><code>airy_bi_diff_scal(z)</code></b>
  </summary>

  - **Description:** Returns a scaled version of the first derivative of the
  Airy function of the second kind and complex argument $z$, i.e.,
  $[\mathrm{d}Bi(z)/\mathrm{d}z] \ e^{-|\mathrm{Re}[(2/3)z^{3/2}]|}$.
  - **Parameter:**
    - `z`, complex argument of
    $[\mathrm{d}Bi(z)/\mathrm{d}z] \ e^{-|\mathrm{Re}[(2/3)z^{3/2}]|}$.
  - **Implementation:** Similar to the `airy_ai()` function.
</details>

## How to use

This library is header‑only, meaning there is nothing to build.

If you still need a compiled file for some reason, refer to the instructions in
the [Compiling the library](#compiling-the-library) section.

Otherwise, you only need to paste all the content of the
[include](include/) folder
inside the include folder of your project (if you do not have an include
folder in your project, paste the content inside the root folder of your
project).

Finally, just write `#include "bessel-library.h"` at the very
beginning of your code and you shall be ready to use the functions.

<details>
  <summary>
    <b>Example of usage in C</b>
  </summary>

```c
#include "bessel-library.h" /* The bessel-library */
#include <complex.h>        /* For double complex type */
#include <stdio.h>          /* For printf() */

int main() {
  double nu = 3.5;
  double complex z = 13.0 + I * 2.7;
  double complex result;
  result = cyl_j(nu, z);
  printf("(%f, %f)\n", creal(result), cimag(result));
  return 0;
}
```
</details>

<details>
  <summary>
    <b>Example of usage in C++</b>
  </summary>

```cpp
#include "bessel-library.h" /* The bessel-library */
#include <complex>          /* For std::complex<double> type */
#include <iostream>         /* For std::cout() */

int main() {
  double nu = 3.5;
  std::complex<double> z = std::complex<double>(13.0, 2.7);
  std::complex<double> result;
  result = cyl_j(nu, z);
  std::cout << "(" << std::real(result) << ", " << std::imag(result) << ")\n";
  return 0;
}
```
</details>

## Some C details

In this library, the implementation is carried out in terms of the C99
standards.

Therefore, all the complex variables are handled using the
`dcomplex` type, which  automatically expands to the `double complex`
type of the C `<complex.h>` library.

Notice that all functions, macros, constants and files whose names contain
the suffix `_impl_` are internal and are not intended to be used by users.

## Compatibility with C++

This library uses `__cplusplus` compiler guards with `extern "C"` and
macros to ensure C++ compatibility (C++98 standard at least).

In this sense, when using C++ compilers, the following C functions are
automatically mapped to their C++ equivalent: `creal(z)`↦`std::real(z)`,
`cimag(z)`↦`std::imag(z)`, `cabs(z)`↦`std::abs(z)`, and
`cexp(z)`↦`std::exp(z)`; and all the complex values are handled by means of
the `dcomplex` type, which automatically expands to the
`std::complex<double>` type of the C++ `<complex>` library, being the macro `I`
the imaginary unity.

## Compiling the library

As aforementioned, usually it is not necessary to compile the library.

However, in any case, the [src](src/) folder contains the file
[bessel-library.c](src/bessel-library.c), which is a C wrapper
([compatible with C++](#compatibility-with-c)) that may be used for
compilation.

The following are examples of how to compile this library using C and C++
compilers.

<details>
  <summary>
    <b>Compiling on Windows with MinGW gcc</b>
  </summary>

  ```bash
  gcc -shared -Iinclude src/bessel-library.c -o bessel-library.dll -Wl,--out-implib,libbessel-library.dll.a
   ```
</details>

<details>
  <summary>
    <b>Compiling on Windows with MinGW g++</b>
  </summary>

  ```bash
  g++ -shared -Iinclude src/bessel-library.c -o bessel-library.dll -Wl,--out-implib,libbessel-library.dll.a
  ```
</details>

<details>
  <summary>
    <b>Compiling on Linux with gcc</b>
  </summary>

  ```bash
  gcc -fPIC -shared -Iinclude src/bessel-library.c -o libbessel-library.so
  ```
</details>

<details>
  <summary>
    <b>Compiling on Linux with g++</b>
  </summary>

  ```bash
  g++ -fPIC -shared -Iinclude src/bessel-library.c -o libbessel-library.so
  ```
</details>

<details>
  <summary>
    <b>Compiling on Windows with MSVC cl</b>
  </summary> 

  Compiling this library with MSVC targeting the C language is discouraged
  because MSVC does not support the `double complex` type from the C
  `<complex.h>` library. However, you may build it in MSVC
  targeting the C++ language (e.g., using the `/TP` flag).

  ```bash
  cl /LD /TP src/bessel-library.c /I include
  ```
</details>

Notice that when using a compiled file of the library into C or C++
projects, you must paste all the content of the
[include](include/) folder inside the include folder of your project,
and then write `#define BESSEL_LIBRARY_IMPORTS`
before the `#include "bessel-library.h"`.

<details>
  <summary>
    <b>Example of usage in C with a compiled file of the library</b>
  </summary>

```c
#define BESSEL_LIBRARY_IMPORTS /* Required for using the compiled file */
#include "bessel-library.h"    /* The bessel-library.h header */
#include <complex.h>           /* For double complex type */
#include <stdio.h>             /* For printf() */

int main() {
  double nu = 3.5;
  double complex z = 13.0 + I * 2.7;
  double complex result;
  result = cyl_j(nu, z);
  printf("(%f, %f)\n", creal(result), cimag(result));
  return 0;
}
```
</details>

## Authorship

The library was developed and is updated by
<a href="https://www.researchgate.net/profile/Jhonas-de-Sarro">
Jhonas O. de Sarro</a> ([@jodesarro](https://github.com/jodesarro)).

Some improvements and tips are also the result of the efforts of the
[contributors](https://github.com/jodesarro/bessel-library/graphs/contributors).

The codes and routines are mainly written based on, or a translation of,
the works listed in the [References](#references) section.

## Licensing

This project is protected under [MIT License](LICENSE).

## References

<a id="ref1"></a>
[1] M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions With
Formulas, Graphs, and Mathematical Tables. Washington, D. C.: National Bureau
of Standards, 1972.

<a id="ref2"></a>
[2] S. Zhang and J. Jin, Computation of Special Functions. New York: Wiley,
1996.

<a id="ref3"></a>
[3] SLATEC Common Mathematical Library, Version 4.1, July 1993. Comprehensive
software library containing over 1400 general purpose mathematical and
statistical routines written in Fortran 77. Available at
<https://www.netlib.org/slatec/> (Accessed: May 25, 2024).