# Bessel library: A C++ library with routines to evaluate Bessel functions of real or complex arguments

<p align="right"><a href="README.pt-br.md">Leia em português (br)</a></p>

## Available features

### <nobr>*J*<sub>0</sub>(*z*)</nobr>
  - Function: `bessel::cyl_j0(z)`.
  - Description: Cylindrical Bessel function of ordinary order for real and complex arguments.
  - Implementation: For real argument, the routine returns the implementation from C++ standard library.
  For complex argument, the routine calculates the ascending series when <nobr>|*z*| ≤ 12</nobr> and Stokes semiconvergent series otherwise with the same idea and truncations presented in Chapter 5 of [[1](#references)] for Fortran.
  Remark that for <nobr>|*z*| > 12</nobr>, `bessel::cyl_j0(z)` strongly depends on the accuracy of functions `std::sin(z)` and `std::cos(z)`.

### New features
  - Unfortunately there are no new features expected anytime soon.

## How to use

The library is in a header-only library style, i.e., there is nothing to build, you only have to include the *<a href="bessel-library.hpp">bessel-library.hpp</a>* file into your project.
See <a href="usage-example.cpp">usage-example.cpp</a> as an example of usage.

## Validity

The function `bessel::cyl_j0(z)` was compared with tabulated values of Chapter 5.9 of [[1](#references)] and with the values of function BesselJ[0,*z*] of [[2](#references)] as one can see in the table below.

|                                      | From [[1](#references)]                                     | From [[2](#references)]                                                    | From `bessel::cyl_j0(z)`                                                          |
|--------------------------------------|-------------------------------------------------------------|----------------------------------------------------------------------------|-----------------------------------------------------------------------------------|
|<nobr>*J*<sub>0</sub>(1+*i*0)</nobr>  |+7.65197687×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>         |+7.651976865579665×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>                 |+7.65197686557966**6**×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>                    |
|<nobr>*J*<sub>0</sub>(5+*i*0)</nobr>  |-1.77596771×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>         |-1.775967713143383×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>                 |-1.77596771314338**4**×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>                    |
|<nobr>*J*<sub>0</sub>(10+*i*0)</nobr> |-2.45935764×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>         |-2.459357644513483×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>                 |-2.459357644513**713**×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>                    |
|<nobr>*J*<sub>0</sub>(25+*i*0)</nobr> |+9.62667833×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>         |+9.626678327595811×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>                 |+9.62667832759581**2**×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>                    |
|<nobr>*J*<sub>0</sub>(50+*i*0)</nobr> |+5.58123277×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>         |+5.581232766925181×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>                 |+5.581232766925181×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>                        |
|<nobr>*J*<sub>0</sub>(100+*i*0)</nobr>|+1.99858503×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>         |+1.998585030422312×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>                 |+1.99858503042231**1**×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>                    |
|<nobr>*J*<sub>0</sub>(4+*i*2)</nobr>  |-1.3787022 ×10<sup>0</sup> <br/>+*i*3.9054236×10<sup>-1</sup>|-1.378702234392483×10<sup>0</sup> <br/>+*i*3.905423570667093×10<sup>-1</sup>|-1.37870223439248**4**×10<sup>0</sup><br/>+*i*3.90542357066709**4**×10<sup>-1</sup>|
|<nobr>*J*<sub>0</sub>(20+*i*10)</nobr>|+1.5460268 ×10<sup>+3</sup><br/>-*i*1.0391216×10<sup>+3</sup>|+1.546026837210333×10<sup>+3</sup><br/>-*i*1.039121575995158×10<sup>+3</sup>|+1.546026837210333×10<sup>+3</sup><br/>-*i*1.039121575995158×10<sup>+3</sup>       |

## Authorship

The codes and routines were developed and are updated by <a href="https://www.researchgate.net/profile/Jhonas-de-Sarro">Jhonas O. de Sarro</a> ([@jodesarro]( https://github.com/jodesarro )).

## Licensing

This project is protected under <a href="LICENSE">MIT License</a>. 

## References

[1] S. Zhang and J. Jin, "Computation of Special Functions", Wiley, 1996.

[2] Wolfram Alpha, accessed 09 April 2022, <https://www.wolframalpha.com/>.
