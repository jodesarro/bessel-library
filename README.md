# Bessel library: A C++ library with routines to evaluate Bessel functions of real or complex arguments

See Refs. [[1–3](#references)] for more information concerning Bessel functions and their computation.

## How to use

The library is in a header-only library style, i.e., there is nothing to build. Therefore, it is very straightforward, you only need to include the <a href="bessel-library.hpp">bessel-library.hpp</a> file into your project (see the usage examples inside the functions listed in the markdown [Available features](#available-features)).

## Available features

The following contains a list of the C++ available functions. Click on each for more information about parameters, implementation, examples, and so on.

### _J_<sub>_ν_</sub>(_z_)—Cylindrical Bessel function of the first kind

<details>
<summary><code>bessel::cyl_j(_nu, _z, _scaled, _flags)</code></summary>
  
- **Description:** Calculation of _cylindrical Bessel function of the first kind_ of integer or real order $\nu$ and real or complex argument $z$, that is, $J_\nu(z)$.
- **Input parameters:**
  - `_nu`: Integer or real order $\nu$ in `T1` type, with `T1` being `int`, `float` or `double`.
  - `_z`: Real argument $z$ in `T2` type, or complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
  - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
  - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
- **Output:** For real $z$, it returns $J_\nu(z)$ in `T2` type; for complex $z$, the complex value of $J_\nu(z)$ in `std::complex<T2>` form.  If `_scaled = true`, it returns $J_\nu(z)\textrm{ }e^{-|\text{Im}(z)|}$.
- **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eqs. (5.4.2) and (5.5.4) of Ref. [[2](#references)] for, respectively, $\nu \in \mathtt{Z}$ and $\nu \notin \mathtt{Z}$; in the latter case, it yields $\infty+i\infty$ when $|z|=0$.
- **Usage example:**
    ```cpp
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
</details>

<details>
  <summary><code>bessel::cyl_j(_nu, _n, _z, _cyl_j, _scaled, _flags)</code></summary>
  
  - **Description:** Concomitant calculation of a number $n \geq 1$ of _cylindrical Bessel functions of the first kind_ of integer or real orders $\nu+k-1$, where $k=1,2,...,n$, and real or complex argument $z$, that is, the sequence $J_\nu(z), J_{\nu+1}(z), ..., J_{\nu+n-1}(z)$.
  - **Input parameters:**
    - `_nu`: Integer or real initial order $\nu$ of the sequence in `T1` type, `T1` being `int`, `float` or `double`.
    - `_n`: Integer greater than zero number $n$ in `int` type.
    - `_z`: Real argument $z$ in `T2` type, or complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
    - `_cyl_j`: Empty array of size `_n`, to sequentially store $J_{\nu+k-1}(z)$ for $k=1,2,...,n$, in `T2` type for real $z$ or in `std::complex<T2>` form for complex $z$.
    - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
    - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
  - **Output:** For real $z$, the real values of $J_\nu(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_j` in `T2` type; for complex $z$, the complex values of $J_{\nu+k-1}(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_j` in `std::complex<T2>` form. If `_scaled = true`, $J_{\nu+k-1}(z)\textrm{ }e^{-|\text{Im}(z)|}$, for $k=1,2,...,n$, are stored.
  - **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eqs. (5.4.2) and (5.5.4) of Ref. [[2](#references)] for, respectively, $\nu \in \mathtt{Z}$ and $\nu \notin \mathtt{Z}$; in the latter case, it yields $\infty+i\infty$ when $|z|=0$.
  - **Usage example:**
    ```cpp
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
</details>

### _Y_<sub>_ν_</sub>(_z_)—Cylindrical Bessel function of the second kind

<details>
  <summary><code>bessel::cyl_y(_nu, _z, _scaled, _flags)</code></summary>
  
  - **Description:** Calculation of _cylindrical Bessel function of the second kind_ of integer or real order $\nu$ and real or complex argument $z$, that is, $Y_\nu(z)$. Such function is also known as _Weber function_ or _Neumann function_, and sometimes written as $N_\nu(z)$.
  - **Input parameters:**
    - `_nu`: Integer or real order $\nu$ in `T1` type, with `T1` being `int`, `float` or `double`.
    - `_z`: Real argument $z$ in `T2` type, or complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
    - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
    - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
  - **Output:** For real $z$, it returns $Y_\nu(z)$ in `T2` type; for complex $z$, the complex value of $Y_\nu(z)$ in `std::complex<T2>` form.  If `_scaled = true`, it returns $Y_\nu(z)\textrm{ }e^{-|\text{Im}(z)|}$.
  - **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eqs. (5.4.2) and (5.5.4) of Ref. [[2](#references)] for, respectively, $\nu \in \mathtt{Z}$ and $\nu \notin \mathtt{Z}$. When $|z|=0$, it yields $-\infty$ if $\nu=0$, or $\infty+i\infty$ otherwise.
  - **Usage example:**
    ```cpp
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
</details>

<details>
  <summary><code>bessel::cyl_y(_nu, _n, _z, _cyl_y, _scaled, _flags)</code></summary>
  
  - **Description:** Concomitant calculation of a number $n \geq 1$ of _cylindrical Bessel functions of the second kind_ of integer or real orders $\nu+k-1$, where $k=1,2,...,n$, and real or complex argument $z$, that is, the sequence $Y_\nu(z), Y_{\nu+1}(z), ..., Y_{\nu+n-1}(z)$. Such functions are also known as _Weber functions_ or _Neumann functions_.
  - **Input parameters:**
    - `_nu`: Integer or real initial order $\nu$ of the sequence in `T1` type, `T1` being `int`, `float` or `double`.
    - `_n`: Integer greater than zero number $n$ in `int` type.
    - `_z`: Real argument $z$ in `T2` type, or complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
    - `_cyl_y`: Empty array of size `_n`, to sequentially store $Y_{\nu+k-1}(z)$ for $k=1,2,...,n$, in `T2` type for real $z$ or in `std::complex<T2>` form for complex $z$.
    - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
    - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
  - **Output:** For real $z$, the real values of $Y_\nu(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_y` in `T2` type; for complex $z$, the complex values of $Y_{\nu+k-1}(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_y` in `std::complex<T2>` form. If `_scaled = true`, $Y_{\nu+k-1}(z)\textrm{ }e^{-|\text{Im}(z)|}$, for $k=1,2,...,n$, are stored.
  - **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eqs. (5.4.2) and (5.5.4) of Ref. [[2](#references)] for, respectively, $\nu \in \mathtt{Z}$ and $\nu \notin \mathtt{Z}$. When $|z|=0$, it yields $-\infty$ if $\nu=0$, or $\infty+i\infty$ otherwise.
  - **Usage example:**
    ```cpp
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
</details>

### _H_<sub>_ν_</sub><sup>(1)</sup>(_z_)—Cylindrical Hankel function of the first kind

<details>
<summary><code>bessel::cyl_h1(_nu, _z, _scaled, _flags)</code></summary>
  
- **Description:** Calculation of _cylindrical Hankel function of the first kind_ of integer or real order $\nu$ and real or complex argument $z$, that is, $H_\nu^{(1)}(z)$. Hankel functions are also known as _Bessel function of the third kind_.
- **Input parameters:**
  - `_nu`: Integer or real order $\nu$ in `T1` type, with `T1` being `int`, `float` or `double`.
  - `_z`: Real argument $z$ in `T2` type, or complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
  - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
  - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
- **Output:** For real $z$, it returns  $H_\nu^{(1)}(z)$ in `T2` type; for complex $z$, the complex value of  $H_\nu^{(1)}(z)$ in `std::complex<T2>` form.  If `_scaled = true`, it returns  $H_\nu^{(1)}(z)\textrm{ }e^{-iz}$.
- **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eq. (9.1.6) of Ref. [[1](#references)]. It yields $\infty+i\infty$ when $|z|=0$.
- **Usage example:**
    ```cpp
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
     result = bessel::cyl_h1( nu, z );
    
     // Alternative calculation with flags
     result = bessel::cyl_h1( nu, z, false, true );
    
     // Calculation of the scaled version
     scaled_result = bessel::cyl_h1( nu, z, true );
    
     // Printing the results
     std::cout << result << std::endl;
     std::cout << scaled_result;
    }
    ```
</details>

<details>
  <summary><code>bessel::cyl_h1(_nu, _n, _z, _cyl_h1, _scaled, _flags)</code></summary>
  
  - **Description:** Concomitant calculation of a number $n \geq 1$ of _cylindrical Hankel functions of the first kind_ of integer or real orders $\nu+k-1$, where $k=1,2,...,n$, and real or complex argument $z$, that is, the sequence $H_\nu^{(1)}(z), H_{\nu+1}^{(1)}(z), ..., H_{\nu+n-1}^{(1)}(z)$. Hankel functions are also known as _Bessel function of the third kind_.
  - **Input parameters:**
    - `_nu`: Integer or real initial order $\nu$ of the sequence in `T1` type, `T1` being `int`, `float` or `double`.
    - `_n`: Integer greater than zero number $n$ in `int` type.
    - `_z`: Real argument $z$ in `T2` type, or complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
    - `_cyl_h1`: Empty array of size `_n`, to sequentially store $H_{\nu+k-1}^{(1)}(z)$ for $k=1,2,...,n$, in `T2` type for real $z$ or in `std::complex<T2>` form for complex $z$.
    - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
    - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
  - **Output:** For real $z$, the real values of $H_{\nu+k-1}^{(1)}(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_h1` in `T2` type; for complex $z$, the complex values of $H_{\nu+k-1}^{(1)}(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_h1` in `std::complex<T2>` form. If `_scaled = true`, $H_{\nu+k-1}^{(1)}(z)\textrm{ }e^{-iz}$, for $k=1,2,...,n$, are stored.
  - **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eq. (9.1.6) of Ref. [[1](#references)]. It yields $\infty+i\infty$ when $|z|=0$.
  - **Usage example:**
    ```cpp
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
      bessel::cyl_h1( nu, 3, z, results );
  
      // Alternative calculation with flags
      bessel::cyl_h1( nu, 3, z, results, false, true );
  
      // Calculation of the scaled versions
      bessel::cyl_h1( nu, 3, z, scaled_results, true );
  
      // Printing the results
      std::cout << results[0] << ", " << results[1] << ", " << results[2] << std::endl;
      std::cout << scaled_results[0] << ", " << scaled_results[1] << ", " << scaled_results[2];
    }
    ```
</details>

### _H_<sub>_ν_</sub><sup>(2)</sup>(_z_)—Cylindrical Hankel function of the second kind

<details>
<summary><code>bessel::cyl_h2(_nu, _z, _scaled, _flags)</code></summary>
  
- **Description:** Calculation of _cylindrical Hankel function of the second kind_ of integer or real order $\nu$ and real or complex argument $z$, that is, $H_\nu^{(2)}(z)$. Hankel functions are also known as _Bessel function of the third kind_.
- **Input parameters:**
  - `_nu`: Integer or real order $\nu$ in `T1` type, with `T1` being `int`, `float` or `double`.
  - `_z`: Real argument $z$ in `T2` type, or complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
  - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
  - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
- **Output:** For real $z$, it returns  $H_\nu^{(2)}(z)$ in `T2` type; for complex $z$, the complex value of  $H_\nu^{(2)}(z)$ in `std::complex<T2>` form.  If `_scaled = true`, it returns  $H_\nu^{(2)}(z)\textrm{ }e^{iz}$.
- **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eq. (9.1.6) of Ref. [[1](#references)]. It yields $\infty+i\infty$ when $|z|=0$.
- **Usage example:**
    ```cpp
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
     result = bessel::cyl_h2( nu, z );
    
     // Alternative calculation with flags
     result = bessel::cyl_h2( nu, z, false, true );
    
     // Calculation of the scaled version
     scaled_result = bessel::cyl_h2( nu, z, true );
    
     // Printing the results
     std::cout << result << std::endl;
     std::cout << scaled_result;
    }
    ```
</details>

<details>
  <summary><code>bessel::cyl_h2(_nu, _n, _z, _cyl_h2, _scaled, _flags)</code></summary>
  
  - **Description:** Concomitant calculation of a number $n \geq 1$ of _cylindrical Hankel functions of the second kind_ of integer or real orders $\nu+k-1$, where $k=1,2,...,n$, and real or complex argument $z$, that is, the sequence $H_\nu^{(2)}(z), H_{\nu+1}^{(2)}(z), ..., H_{\nu+n-1}^{(2)}(z)$. Hankel functions are also known as _Bessel function of the third kind_.
  - **Input parameters:**
    - `_nu`: Integer or real initial order $\nu$ of the sequence in `T1` type, `T1` being `int`, `float` or `double`.
    - `_n`: Integer greater than zero number $n$ in `int` type.
    - `_z`: Real argument $z$ in `T2` type, or complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
    - `_cyl_h2`: Empty array of size `_n`, to sequentially store $H_{\nu+k-1}^{(2)}(z)$ for $k=1,2,...,n$, in `T2` type for real $z$ or in `std::complex<T2>` form for complex $z$.
    - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
    - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
  - **Output:** For real $z$, the real values of $H_{\nu+k-1}^{(2)}(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_h2` in `T2` type; for complex $z$, the complex values of $H_{\nu+k-1}^{(2)}(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_h2` in `std::complex<T2>` form. If `_scaled = true`, $H_{\nu+k-1}^{(2)}(z)\textrm{ }e^{iz}$, for $k=1,2,...,n$, are stored.
  - **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eq. (9.1.6) of Ref. [[1](#references)]. It yields $\infty+i\infty$ when $|z|=0$.
  - **Usage example:**
    ```cpp
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
      bessel::cyl_h2( nu, 3, z, results );
  
      // Alternative calculation with flags
      bessel::cyl_h2( nu, 3, z, results, false, true );
  
      // Calculation of the scaled versions
      bessel::cyl_h2( nu, 3, z, scaled_results, true );
  
      // Printing the results
      std::cout << results[0] << ", " << results[1] << ", " << results[2] << std::endl;
      std::cout << scaled_results[0] << ", " << scaled_results[1] << ", " << scaled_results[2];
    }
    ```
</details>

### _I_<sub>_ν_</sub>(_z_)—Modified cylindrical Bessel function of the first kind

<details>
<summary><code>bessel::cyl_i(_nu, _z, _scaled, _flags)</code></summary>
  
- **Description:** Calculation of _modified cylindrical Bessel function of the first kind_ of integer or real order $\nu$ and real or complex argument $z$, that is, $I_\nu(z)$. Such function is also known as _cylindrical Bessel function of imaginary argument_ or sometimes as _hyperbolic Bessel function_.
- **Input parameters:**
  - `_nu`: Integer or real order $\nu$ in `T1` type, with `T1` being `int`, `float` or `double`.
  - `_z`: Real argument $z$ in `T2` type, or complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
  - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
  - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
- **Output:** For real $z$, it returns $I_\nu(z)$ in `T2` type; for complex $z$, the complex value of $I_\nu(z)$ in `std::complex<T2>` form.  If `_scaled = true`, it returns $I_\nu(z)\textrm{ }e^{-|\text{Re}(z)|}$.
- **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eqs. (6.1.5) and (6.5.4) of Ref. [[2](#references)] for, respectively, $\nu \in \mathtt{Z}$ and $\nu \notin \mathtt{Z}$; in the latter case, it yields $\infty+i\infty$ when $|z|=0$.
- **Usage example:**
    ```cpp
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
     result = bessel::cyl_i( nu, z );
    
     // Alternative calculation with flags
     result = bessel::cyl_i( nu, z, false, true );
    
     // Calculation of the scaled version
     scaled_result = bessel::cyl_i( nu, z, true );
    
     // Printing the results
     std::cout << result << std::endl;
     std::cout << scaled_result;
    }
    ```
</details>

<details>
  <summary><code>bessel::cyl_i(_nu, _n, _z, _cyl_i, _scaled, _flags)</code></summary>
  
  - **Description:** Concomitant calculation of a number $n \geq 1$ of _modified cylindrical Bessel functions of the first kind_ of integer or real orders $\nu+k-1$, where $k=1,2,...,n$, and real or complex argument $z$, that is, the sequence $I_\nu(z), I_{\nu+1}(z), ..., I_{\nu+n-1}(z)$. Such functions are also known as _cylindrical Bessel functions of imaginary argument_ or sometimes as _hyperbolic Bessel functions_.
  - **Input parameters:**
    - `_nu`: Integer or real initial order $\nu$ of the sequence in `T1` type, `T1` being `int`, `float` or `double`.
    - `_n`: Integer greater than zero number $n$ in `int` type.
    - `_z`: Real argument $z$ in `T2` type, or complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
    - `_cyl_i`: Empty array of size `_n`, to sequentially store $I_{\nu+k-1}(z)$ for $k=1,2,...,n$, in `T2` type for real $z$ or in `std::complex<T2>` form for complex $z$.
    - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
    - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
  - **Output:** For real $z$, the real values of $I_\nu(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_i` in `T2` type; for complex $z$, the complex values of $I_{\nu+k-1}(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_i` in `std::complex<T2>` form. If `_scaled = true`, $I_{\nu+k-1}(z)\textrm{ }e^{-|\text{Re}(z)|}$, for $k=1,2,...,n$, are stored.
  - **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eqs. (6.1.5) and (6.5.4) of Ref. [[2](#references)] for, respectively, $\nu \in \mathtt{Z}$ and $\nu \notin \mathtt{Z}$; in the latter case, it yields $\infty+i\infty$ when $|z|=0$.
  - **Usage example:**
    ```cpp
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
      bessel::cyl_i( nu, 3, z, results );
  
      // Alternative calculation with flags
      bessel::cyl_i( nu, 3, z, results, false, true );
  
      // Calculation of the scaled versions
      bessel::cyl_i( nu, 3, z, scaled_results, true );
  
      // Printing the results
      std::cout << results[0] << ", " << results[1] << ", " << results[2] << std::endl;
      std::cout << scaled_results[0] << ", " << scaled_results[1] << ", " << scaled_results[2];
    }
    ```
</details>

### _K_<sub>_ν_</sub>(_z_)—Modified cylindrical Bessel function of the second kind

<details>
  <summary><code>bessel::cyl_k(_nu, _z, _scaled, _flags)</code></summary>
  
  - **Description:** Calculation of _modified cylindrical Bessel function of the second kind_ of integer or real order $\nu$ and real or complex argument $z$, that is, $K_\nu(z)$. Such function is also known as _Basset function_ or _MacDonald function_.
  - **Input parameters:**
    - `_nu`: Integer or real order $\nu$ in `T1` type, with `T1` being `int`, `float` or `double`.
    - `_z`: Real argument $z$ in `T2` type, or complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
    - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
    - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
  - **Output:** For real $z$, it returns $K_\nu(z)$ in `T2` type; for complex $z$, the complex value of $K_\nu(z)$ in `std::complex<T2>` form.  If `_scaled = true`, it returns $K_\nu(z)\textrm{ }e^{z}$.
  - **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eqs. (6.5.5) of Ref. [[2](#references)]. When $|z|=0$, it yields $\infty$ if $\nu=0$, or $\infty+i\infty$ otherwise.
  - **Usage example:**
    ```cpp
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
      result = bessel::cyl_k( nu, z );
    
      // Alternative calculation with flags
      result = bessel::cyl_k( nu, z, false, true );
  
      // Calculation of the scaled version
      scaled_result = bessel::cyl_k( nu, z, true );
      
      // Printing the results
      std::cout << result << std::endl;
      std::cout << scaled_result;
    }
    ```
</details>

<details>
  <summary><code>bessel::cyl_k(_nu, _n, _z, _cyl_k, _scaled, _flags)</code></summary>

  - **Description:** Concomitant calculation of a number $n \geq 1$ of _modified cylindrical Bessel functions of the second kind_ of integer or real orders $\nu+k-1$, where $k=1,2,...,n$, and real or complex argument $z$, that is, the sequence $K_\nu(z), K_{\nu+1}(z), ..., K_{\nu+n-1}(z)$. Such functions are also known as _Basset functions_ or _MacDonald functions_.
  - **Input parameters:**
    - `_nu`: Integer or real initial order $\nu$ of the sequence in `T1` type, `T1` being `int`, `float` or `double`.
    - `_n`: Integer greater than zero number $n$ in `int` type.
    - `_z`: Real argument $z$ in `T2` type, or complex argument $z$ in `std::complex<T2>` form, with `T2` being `float` or `double`.
    - `_cyl_k`: Empty array of size `_n`, to sequentially store $K_{\nu+k-1}(z)$ for $k=1,2,...,n$, in `T2` type for real $z$ or in `std::complex<T2>` form for complex $z$.
    - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
    - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
  - **Output:** For real $z$, the real values of $K_\nu(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_k` in `T2` type; for complex $z$, the complex values of $K_{\nu+k-1}(z)$, for $k=1,2,...,n$, are stored in the array `_cyl_k` in `std::complex<T2>` form. If `_scaled = true`, $K_{\nu+k-1}(z)\textrm{ }e^{z}$, for $k=1,2,...,n$, are stored.
  - **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library. Negative orders are handled by Eqs. (6.5.5) of Ref. [[2](#references)]. When $|z|=0$, it yields $\infty$ if $\nu=0$, or $\infty+i\infty$ otherwise.
  - **Usage example:**
    ```cpp
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
      bessel::cyl_k( nu, 3, z, results );

      // Alternative calculation with flags
      bessel::cyl_k( nu, 3, z, results, false, true );

      // Calculation of the scaled versions
      bessel::cyl_k( nu, 3, z, scaled_results, true );

      // Printing the results
      std::cout << results[0] << ", " << results[1] << ", " << results[2] << std::endl;
      std::cout << scaled_results[0] << ", " << scaled_results[1] << ", " << scaled_results[2];
    }
    ```
</details>

### Ai(_z_)—Airy function of the first kind

<details>
  <summary><code>bessel::airy_ai(_z, _scaled, _flags)</code></summary>
  
  - **Description:** Calculation of _Airy function of the first kind_ of real or complex argument $z$, that is, $\textrm{Ai}(z)$.
  - **Input parameters:**
    - `_z`: Real argument $z$ in `T1` type, or complex argument $z$ in `std::complex<T1>` form, with `T1` being `float` or `double`.
    - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
    - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
  - **Output:** For real $z$, it returns $\textrm{Ai}(z)$ in `T1` type; for complex $z$, the complex value of $\textrm{Ai}(z)$ in `std::complex<T1>` form.  If `_scaled = true`, it returns $\textrm{Ai}(z)\textrm{ }e^{(2/3)z\sqrt{z}}$.
  - **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library.
  - **Usage example:**
    ```cpp
    #include <iostream>
    #include "bessel-library.hpp"
  
    int main()
    {
      // Declaration of variables
      std::complex<double> z = std::complex<double>(1.0, 2.3);
      std::complex<double> result;
      std::complex<double> scaled_result;
    
      // Calculation of a function
      result = bessel::airy_ai( z );
    
      // Alternative calculation with flags
      result = bessel::airy_ai( z, false, true );
  
      // Calculation of the scaled version
      scaled_result = bessel::airy_ai( z, true );
      
      // Printing the results
      std::cout << result << std::endl;
      std::cout << scaled_result;
    }
    ```
</details>

### Bi(_z_)—Airy function of the second kind

<details>
  <summary><code>bessel::airy_bi(_z, _scaled, _flags)</code></summary>
  
  - **Description:** Calculation of _Airy function of the second kind_ of real or complex argument $z$, that is, $\textrm{Bi}(z)$.
  - **Input parameters:**
    - `_z`: Real argument $z$ in `T1` type, or complex argument $z$ in `std::complex<T1>` form, with `T1` being `float` or `double`.
    - `_scaled`: Optional `bool` parameter. If `true`, returns a scaled version of the result (see Output topic below).
    - `_flags`: Optional `bool` parameter. If `true`, print error and warning messages.
  - **Output:** For real $z$, it returns $\textrm{Bi}(z)$ in `T1` type; for complex $z$, the complex value of $\textrm{Bi}(z)$ in `std::complex<T1>` form.  If `_scaled = true`, it returns $\textrm{Bi}(z)\textrm{ }e^{-|\textrm{Re}[(2/3)z\sqrt{z}]|}$.
  - **Implementation:** In general, the routine is based on the D. E. Amos Fortran 77 routines of the SLATEC library [[3](#references)]. Such Fortran routines, and all their dependencies, were carefully translated to be used in this library.
  - **Usage example:**
    ```cpp
    #include <iostream>
    #include "bessel-library.hpp"
  
    int main()
    {
      // Declaration of variables
      std::complex<double> z = std::complex<double>(1.0, 2.3);
      std::complex<double> result;
      std::complex<double> scaled_result;
    
      // Calculation of a function
      result = bessel::airy_bi( z );
    
      // Alternative calculation with flags
      result = bessel::airy_bi( z, false, true );
  
      // Calculation of the scaled version
      scaled_result = bessel::airy_bi( z, true );
      
      // Printing the results
      std::cout << result << std::endl;
      std::cout << scaled_result;
    }
    ```
</details>

## New features coming soon
  - Implementation of spherical versions of the functions.

## Authorship

The codes and routines were developed and are updated by <a href="https://www.researchgate.net/profile/Jhonas-de-Sarro">Jhonas O. de Sarro</a> ([@jodesarro]( https://github.com/jodesarro )). They are mainly written based on, or a translation of, the content of Refs. [[1–3](#references)].

## Licensing

This project is protected under <a href="LICENSE">MIT License</a>. 

## References

[1] M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions With Formulas, Graphs, and Mathematical Tables. Washington, D. C.: National Bureau of Standards, 1972.<br/>
[2] S. Zhang and J. Jin, Computation of Special Functions. New York: Wiley, 1996.<br/>
[3] SLATEC Common Mathematical Library, Version 4.1, July 1993. Comprehensive software library containing over 1400 general purpose mathematical and statistical routines written in Fortran 77. Available at https://www.netlib.org/slatec/ (Accessed: May 25, 2024).<br/>