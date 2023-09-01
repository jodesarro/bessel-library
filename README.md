# Bessel library: A C++ library with routines to evaluate Bessel functions of real or complex arguments

<p align="right"><a href="README.pt-br.md">Leia em português (br)</a></p>

## Available features

### <nobr>*J*<sub>n</sub>(*z*)</nobr>
  - **Function:** `bessel::cyl_j(n,z)`.
  - **Description:** Cylindrical Bessel function of the first kind for an integer order *n* and a real or complex argument *z*.
  - **Implementation for real arguments:** The routine returns the implementation from C++ standard library.
  - **Implementation for complex arguments:** For <nobr>|*n*| ≤ 1</nobr>, the routine calculates the ascending series if <nobr>|*z*| ≤ 12</nobr> or the Stokes semiconvergent series otherwise.
  For <nobr>|*n*| > 1</nobr>, the routine makes use of a forward recurrence if <nobr>*n* < |*z*|/4</nobr> or a backward recurrence otherwise.
  The implementation follows the same idea as the routines presented in Chapter 5 of <nobr>Ref. [[1](#references)]</nobr> for Fortran.

### <nobr>*J*<sub>n</sub>'(*z*)</nobr>
  - **Function:** `bessel::cyl_j_diff(n,z)`.
  - **Description:** Derivative of <nobr>*J*<sub>n</sub>(*z*)</nobr> with respect to *z*.
  - **Implementation for real arguments:** The routine uses a recurrence relation together with the implementation from C++ standard library.
  - **Implementation for complex arguments:** The routine returns `-bessel::cyl_j(1,z)` if <nobr>*n* = 0</nobr> or uses `bessel::cyl_j(n,z)` in a recurrence relation otherwise.

### New features
  - Unfortunately there are no new features expected anytime soon.

## How to use

The library is in a header-only library style, i.e., there is nothing to build, you only have to include the *<a href="bessel-library.hpp">bessel-library.hpp</a>* file into your project.
See <a href="usage-example.cpp">usage-example.cpp</a> as an example of usage.

## Validity

For some given parameters, the functions were compared with the function BesselJ[*n*,*z*] of <nobr>Ref. [[2](#references)]</nobr>. The results are listed below.

|                                            | From Ref. [[2](#references)]                                                   | From this library                                                                    |
|--------------------------------------------|--------------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
|<nobr>*J*<sub>0</sub>(3.2-*i*1.7)</nobr>    |-1.030282834527403×10<sup>0</sup><br/>  +*i*5.419144653871709×10<sup>-1</sup>   |-1.03028283452740**4**×10<sup>0</sup><br/>  +*i*5.41914465387170**5**×10<sup>-1</sup> |
|<nobr>*J*<sub>0</sub>(-15.1-*i*83.9)</nobr> |-9.079800270662861×10<sup>+34</sup><br/>-*i*7.605394317402608×10<sup>+34</sup>  |-9.079800270662**912**×10<sup>+34</sup><br/>-*i*7.6053943174026**56**×10<sup>+34</sup>|
|<nobr>*J*<sub>-1</sub>(-2.8+*i*7.7)</nobr>  |+1.405657382170826×10<sup>+2</sup><br/> +*i*2.583924513755478×10<sup>+2</sup>   |+1.40565738217082**8**×10<sup>+2</sup><br/> +*i*2.58392451375547**9**×10<sup>+2</sup> |
|<nobr>*J*<sub>3</sub>(4+*i*29)</nobr>       |+1.807884960557100×10<sup>+11</sup><br/>+*i*1.717728203190626×10<sup>+11</sup>  |+1.807884960557100×10<sup>+11</sup><br/>    +*i*1.717728203190626×10<sup>+11</sup>    |
|<nobr>*J*<sub>-14</sub>(42-*i*8)</nobr>     |-1.109381836657852×10<sup>+2</sup><br/> +*i*4.649254392531784×10<sup>+1</sup>   |-1.109381836657852×10<sup>+2</sup><br/>     +*i*4.64925439253178**5**×10<sup>+1</sup> |
|<nobr>*J*<sub>100</sub>(71+*i*30)</nobr>    |-1.873523256314857×10<sup>-4</sup><br/> +*i*2.688597190207697×10<sup>-6</sup>   |-1.8735232563148**47**×10<sup>-4</sup><br/> +*i*2.68859719020**8606**×10<sup>-6</sup> |
|<nobr>*J*<sub>0</sub>'(-11+*i*18)</nobr>    |-5.427270167814426×10<sup>+6</sup><br/> +*i*1.444852021146701×10<sup>+6</sup>   |-5.427270167814426×10<sup>+6</sup><br/>     +*i*1.44485202114670**2**×10<sup>+6</sup> |
|<nobr>*J*<sub>-9</sub>'(-6.6-*i*3.6)</nobr> |+1.608080761747265×10<sup>-1</sup><br/> -*i*1.505778965628776×10<sup>-1</sup>   |+1.60808076174726**4**×10<sup>-1</sup><br/> -*i*1.50577896562877**1**×10<sup>-1</sup> |
|<nobr>*J*<sub>47</sub>'(2.71+*i*3.14)</nobr>|-4.836547831119548×10<sup>-45</sup><br/>+*i*3.403725372134759×10<sup>-44</sup>  |-4.8365478311195**86**×10<sup>-45</sup><br/>+*i*3.40372537213475**6**×10<sup>-44</sup>|

## Authorship

The codes and routines were developed and are updated by <a href="https://www.researchgate.net/profile/Jhonas-de-Sarro">Jhonas O. de Sarro</a> ([@jodesarro]( https://github.com/jodesarro )).

## Licensing

This project is protected under <a href="LICENSE">MIT License</a>. 

## References

[1] S. Zhang and J. Jin, "Computation of Special Functions", Wiley, 1996.

[2] Wolfram Alpha, accessed 01 September 2023, <https://www.wolframalpha.com/>.