# Bessel library: A C++ library with routines to evaluate Bessel functions of real or complex arguments

<p align="right"><a href="README.pt-br.md">Leia em português (br)</a></p>

## Available features

### <nobr>*J*<sub>n</sub>(*z*)</nobr>
  - **Function:** `bessel::cyl_j(n,z)`.
  - **Description:** Cylindrical Bessel function of the first kind for an integer order *n* and a real or complex argument *z*.
  - **Options:** The function has a third optional boolean argument called *warnings*: `bessel::cyl_j(n,z,warnings)`. The standard value for *warnings* is `true`, set it to `false` to disable warnings. 
  - **Implementation for real arguments:** The routine returns the implementation from C++ standard library.
  - **Implementation for complex arguments:** The routine is based in the calculation of ascending series, Stokes semiconvergent series, forward recurrences and backward recurrences.
  The implementation follows the same idea as the routines presented in Chapter 5 of <nobr>Ref. [[1](#references)]</nobr> for Fortran.

### <nobr>*J*<sub>n</sub>'(*z*)</nobr>
  - **Function:** `bessel::cyl_j_diff(n,z)`.
  - **Description:** Derivative of <nobr>*J*<sub>n</sub>(*z*)</nobr> with respect to *z*.
  - **Options:** The function has a third optional boolean argument called *warnings*: `bessel::cyl_j_diff(n,z,warnings)`. The standard value for *warnings* is `true`, set it to `false` to disable warnings. 
  - **Implementation for real arguments:** The routine uses a recurrence relation together with the implementation from C++ standard library.
  - **Implementation for complex arguments:** The routine returns `-bessel::cyl_j(1,z)` if <nobr>*n* = 0</nobr> or uses `bessel::cyl_j(n,z)` in a recurrence relation otherwise.

### New features
  - Unfortunately there are no new features expected anytime soon.

## How to use

The library is in a header-only library style, i.e., there is nothing to build, you only have to include the *<a href="bessel-library.hpp">bessel-library.hpp</a>* file into your project.
See <a href="usage-example.cpp">usage-example.cpp</a> as an example of usage.

## Validity

The routines were written to achieve results with single precision for <nobr>|Im(*z*)| ≤ 21</nobr> when compared with the function BesselJ[*n*,*z*] of <nobr>Ref. [[2](#references)]</nobr>.

## Authorship

The codes and routines were developed and are updated by <a href="https://www.researchgate.net/profile/Jhonas-de-Sarro">Jhonas O. de Sarro</a> ([@jodesarro]( https://github.com/jodesarro )).

## Licensing

This project is protected under <a href="LICENSE">MIT License</a>. 

## References

[1] S. Zhang and J. Jin, "Computation of Special Functions", Wiley, 1996.

[2] Wolfram Alpha, accessed 01 September 2023, <https://www.wolframalpha.com/>.