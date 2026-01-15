# CHANGELOG

Refer to the reference section of the [README.md](README.md) file
for a complete list of references.

## v1.0.1 Jan 13, 2026

- In Airy functions, all `_deriv` terms were changed to `_diff`.
- Include paths corrected in file include/bessel-library.h
- For increments standardized as i++.
- More descriptions were included in the comments for some functions.
- Added macro for `static inline` for all core functions. Now
src/bessel-library.c just overwrites such macro.
- Added file src/bessel-library-declarations.c containing a list of all the
functions declarations.
- Inclusion of the header version.h in the include/bessel-library.h file.
- A copy of the license was included in include/bessel-library/license.txt,
and reference to it were included at the beginning of all code files.
- README.me updated.

## v1.0.0 Jan 09, 2026

- The library was consistently renamed to "Bessel Library: A C library with
routines for computing Bessel functions".
- The library was rewritten in C using the C99 standards, but is still
fully compatible with C++ (C++98 standards) through the use of "__cplusplus"
guards. See the sections "Some C details" and "Compatibility with C++" of the
README.md file.
- The library was splitted in several files and folders. The folder impl
contains internal routines not intended to be called by the users. All 
functions, macros, constants and files whose names contain the suffix `_impl_`
are internal and are not intended to be used by users. The folder core
contains isolated headers for functions, and they may be called by users.
- Some functions and their parameters were consistently adequated to C, as
described in the section "Available features" of the README.md file.
- Scaled version of I_nu(z) for nu<0 was corrected based on K_nu(z).
- README.md file was consistently rewritten.
- Creation of the version.h file for version control.
- Creation of the references.txt file for a list of the references.
- Creation of the CHANGELOG.md file for version control.
- Standardization of all beginnings of .h files.
- Creation a bessel-library.c file inside the src folder for compilation
purposes.

## v0.1.3 Jun 09, 2024

- Inclusion of Hankel functions cyl_h1 and cyl_h2 for integer or real orders
and real or complex arguments.
- Inclusion of Airy functions airy_ai and airy_bi for real or complex
arguments.
- Functions mod_i and mod_k were consistently renamed to cyl_i and cyl_k.
- The flags now print the number of components set to zero due to underflow.
- Routines zairy_, zbesh_, zbesj_, zbesy_, zbesi_, zbesk_ and zbiry_, based
on version 930101 of D. E. Amos routines (https://doi.org/10.1145/7921.214331,
ACM domain) were changed (reverted) to be based on slatec ([3], public domain)
versions to avoid copyright conflicts between ACM and MIT licenses and
permissions. The versions v0.1.1 and v0.1.2 of this code, and github commits
related to them, shall be deleted and must be disconsidered and discarded by
all users.
- Revision and reorganization of all slatec functions.
- Creation of functions d1mach and i1mach to make easier to compare with
original slatec versions.

## v0.1.2 Jun 06, 2024

- Inclusion of modified Bessel functions mod_i and mod_k for integer or real
orders and real or complex arguments.

## v0.1.1 May 27, 2024

- Routines zairy_, zbesh_, zbesj_, and zbesy_, updated to the version 930101
of D. E. Amos routines (https://doi.org/10.1145/7921.214331).
- Inclusion of routines zbesi_, zbesk_, zbiry_ accordingly to version 930101
of D. E. Amos routines (https://doi.org/10.1145/7921.214331).
- Inclusion of C++ callable functions to overload cyl_j and cyl_y for real
arguments.
- Static declarations removed for thread safety.

## v0.1.0 May 26, 2024

- Routines for cyl_j based on Ref. [2] were replaced by D. E. Amos Fortran 77
routines of SLATEC library [3].
- D. E. Amos routines zairy_.f, zbesh_.f, zbesj_.f, zbesy_.f, and all their
dependencies, were converted to C using f2c (Availabe at:
https://www.netlib.org/f2c/. Accessed: May 25, 2024).
- Replacement of all functions d1mach amd i1mach by C macros of float.h.
- Corrections of the translated f2c version and elimination of external
dependencies.
- Reorganization of the whole code to be easily callable from C++.
- Inclusion of cylindrical Bessel functions of the second kind
(or Neumann functions) cyl_y.
- Calculation of negative orders for cyl_j and cyl_y through Eq. (5.5.4) of
Ref. [2].
- Now, cyl, Bessel functions of the first and second kinds, cyl_j and cyl_y,
are available also for real (positive or negative) orders.
- Inclusion of cyl_j and cyl_y that returns an array of an int sequence of
orders.
- Inclusion of parameters to print flag messages, and to return scaled
versions of cyl_j and cyl_y.
- Inclusion of namespace bessel::slatec to call all slatec routines.

## v0.0.0 until May 12, 2024 

- Routines for cylindrical Bessel functions of the first kind and int order
written based on Ref. [2].