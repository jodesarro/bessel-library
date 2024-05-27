# Bessel library: A C++ library with routines to evaluate Bessel functions of real or complex arguments

See Refs. [[1–3](#references)] for more information concerning Bessel functions.

## Available features

### $J_\nu(z)$—Cylindrical Bessel function of the first kind

#### Function `bessel::cyl_j(_nu, _z, _scaled, _flags)`
- **Description:** Calculation of cylindrical Bessel function of the first kind of real order $\nu$ and complex argument $z$, that is, $J_\nu(z)$.
- **Input parameters:**
  - `_nu`: real order $\nu$ in `T1` type, `T1` being `float` or `double`.
  - `_z`: complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
  - `_scaled`: optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
  - `_flags`: optional `bool` parameter. If `true`, print error and warning messages.
- **Output:** It returns the complex value of $J_\nu(z)$ in `std::complex<T2>` form; if `_scaled = true`, it returns $J_\nu(z)\textrm{ }e^{-|\text{Im}(z)|}$.
- **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[2](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eq. (5.5.4) of Ref. [[1](#references)]. For $|z|=0$ with $\nu \notin \mathrm{Z}$, it gives $\infty+i\infty$.
- **Usage example:**
  ```
  #include <iostream>
  #include "bessel-library.hpp"

  int main()
  {
    // Declaration of variables
    double nu = 3.5;
    std::complex<double> z = std::complex<double>(1.0, 2.3);
    std::complex<double> result;
    std::complex<double> scaled_result;
  
    // Calculation of a function
    result = bessel::cyl_j( nu, z );
  
    // Alternative calculation with flags
    result = bessel::cyl_j( nu, z, false, true );

    // Calculation of the scaled version
    scaled_result = bessel::cyl_j( nu, z, true );
    
    // Printing the results
    std::cout << result << std::endl;
    std::cout << scaled_result;
  }
  ```

#### Function `bessel::cyl_j(_nu, _n, _z, _cyl_j, _scaled, _flags)`
- **Description:** Concomitant calculation of a number $n$ of cylindrical Bessel functions of the first kind of real orders $\nu+k-1$, where $k=1,2,...,n$, and complex argument $z$, that is, the sequence $J_\nu(z), J_{\nu+1}(z), ..., J_{\nu+n-1}(z)$.
- **Input parameters:**
  - `_nu`: real initial order $\nu$ of the sequence in `T1` type, `T1` being `float` or `double`.
  - `_n`: integer number $n$ in `int` type.
  - `_z`: complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
  - `_cyl_j`: empty array of size `_n` to sequentially store $J_{\nu+k-1}(z)$ for $k=1,2,...,n$ in `std::complex<T2>` form.
  - `_scaled`: optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
  - `_flags`: optional `bool` parameter. If `true`, print error and warning messages.
- **Output:** The complex values of $J_{\nu+k-1}(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_j` in `std::complex<T2>` form; if `_scaled = true`, $J_{\nu+k-1}(z)\textrm{ }e^{-|\text{Im}(z)|}$, for $k=1,2,...,n$, are stored.
- **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[2](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eq. (5.5.4) of Ref. [[1](#references)]. For $|z|=0$ with $\nu \notin \mathrm{Z}$, it gives $\infty+i\infty$.
- **Usage example:**
  ```
  #include <iostream>
  #include "bessel-library.hpp"

  int main()
  {
    // Declaration of variables
    double nu = -1.7;
    std::complex<double> z = std::complex<double>(1.2, 5.3);
    std::complex<double> results [3];
    std::complex<double> scaled_results [3];
    
    // Calculation of functions
    bessel::cyl_j( nu, 3, z, results );

    // Alternative calculation with flags
    bessel::cyl_j( nu, 3, z, results, false, true );

    // Calculation of the scaled versions
    bessel::cyl_j( nu, 3, z, scaled_results, true );

    // Printing the results
    std::cout << results[0] << ", " << results[1] << ", " << results[2] << std::endl;
    std::cout << scaled_results[0] << ", " << scaled_results[1] << ", " << scaled_results[2];
  }
  ```

### $Y_\nu(z)$—Cylindrical Bessel function of the second kind

#### Function `bessel::cyl_y(_nu, _z, _scaled, _flags)`
- **Description:** Calculation of cylindrical Bessel function of the second kind (also known as Neumann function) of real order $\nu$ and complex argument $z$, that is, $Y_\nu(z)$ [also written as $N_\nu(z)$].
- **Input parameters:**
  - `_nu`: real order $\nu$ in `T1` type, `T1` being `float` or `double`.
  - `_z`: complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
  - `_scaled`: optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
  - `_flags`: optional `bool` parameter. If `true`, print error and warning messages.
- **Output:** It returns the complex value of $Y_\nu(z)$ in `std::complex<T2>` form; if `_scaled = true`, it returns $Y_\nu(z)\textrm{ }e^{-|\text{Im}(z)|}$.
- **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[2](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eq. (5.5.4) of Ref. [[1](#references)]. For $|z|=0$ with $\nu \notin \mathrm{Z}$, it gives $\infty+i\infty$.
- **Usage example:**
  ```
  #include <iostream>
  #include "bessel-library.hpp"

  int main()
  {
    // Declaration of variables
    double nu = 3.5;
    std::complex<double> z = std::complex<double>(1.0, 2.3);
    std::complex<double> result;
    std::complex<double> scaled_result;
  
    // Calculation of a function
    result = bessel::cyl_y( nu, z );
  
    // Alternative calculation with flags
    result = bessel::cyl_y( nu, z, false, true );

    // Calculation of the scaled version
    scaled_result = bessel::cyl_y( nu, z, true );
    
    // Printing the results
    std::cout << result << std::endl;
    std::cout << scaled_result;
  }
  ```

#### Function `bessel::cyl_y(_nu, _n, _z, _cyl_y, _scaled, _flags)`
- **Description:** Concomitant calculation of a number $n$ of cylindrical Bessel functions of the second kind (also known as Neumann functions) of real orders $\nu+k-1$, where $k=1,2,...,n$, and complex argument $z$, that is, the sequence $Y_\nu(z), Y_{\nu+1}(z), ..., Y_{\nu+n-1}(z)$.
- **Input parameters:**
  - `_nu`: real initial order $\nu$ of the sequence in `T1` type, `T1` being `float` or `double`.
  - `_n`: integer number $n$ in `int` type.
  - `_z`: complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
  - `_cyl_y`: empty array of size `_n` to sequentially store $Y_{\nu+k-1}(z)$ for $k=1,2,...,n$ in `std::complex<T2>` form.
  - `_scaled`: optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
  - `_flags`: optional `bool` parameter. If `true`, print error and warning messages.
- **Output:** The complex values of $Y_{\nu+k-1}(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_y` in `std::complex<T2>` form; if `_scaled = true`, $Y_{\nu+k-1}(z)\textrm{ }e^{-|\text{Im}(z)|}$, for $k=1,2,...,n$, are stored.
- **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[2](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eq. (5.5.4) of Ref. [[1](#references)]. For $|z|=0$ with $\nu \notin \mathrm{Z}$, it gives $\infty+i\infty$.
- **Usage example:**
  ```
  #include <iostream>
  #include "bessel-library.hpp"

  int main()
  {
    // Declaration of variables
    double nu = -1.7;
    std::complex<double> z = std::complex<double>(1.2, 5.3);
    std::complex<double> results [3];
    std::complex<double> scaled_results [3];
    
    // Calculation of functions
    bessel::cyl_y( nu, 3, z, results );

    // Alternative calculation with flags
    bessel::cyl_y( nu, 3, z, results, false, true );

    // Calculation of the scaled versions
    bessel::cyl_y( nu, 3, z, scaled_results, true );

    // Printing the results
    std::cout << results[0] << ", " << results[1] << ", " << results[2] << std::endl;
    std::cout << scaled_results[0] << ", " << scaled_results[1] << ", " << scaled_results[2];
  }
  ```

### New features coming soon
  - Implementation of Hankel functions.
  - Implementation of Airy functions.

## How to use

The library is in a header-only library style, i.e., there is nothing to build, you only have to include the *<a href="bessel-library.hpp">bessel-library.hpp</a>* file into your project.

## Authorship

The codes and routines were developed and are updated by <a href="https://www.researchgate.net/profile/Jhonas-de-Sarro">Jhonas O. de Sarro</a> ([@jodesarro]( https://github.com/jodesarro )).

## Licensing

This project is protected under <a href="LICENSE">MIT License</a>. 

## References

[1] M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions With Formulas, Graphs, and Mathematical Tables. Washington, D. C.: National Bureau of Standards, 1972.<br/>
[2] S. Zhang and J. Jin, Computation of Special Functions. New York: Wiley, 1996.<br/>
[3] SLATEC Common Mathematical Library, Version 4.1, July 1993. Comprehensive software library containing over 1400 general purpose mathematical and statistical routines written in Fortran 77. Available at https://www.netlib.org/slatec/ (Accessed: May 25, 2024).