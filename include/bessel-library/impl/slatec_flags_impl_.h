/* 
    Bessel Library: A C library with routines for computing Bessel functions

    File: include/bessel-library/impl/slatec_flags_impl_.h
    Version: include/bessel-library/version.h
    Author: Jhonas Olivati de Sarro
    Language standards: C99
    References: include/bessel-library/references.txt
    License: include/bessel-library/license.txt

    Description:
        Functions for printing in stderr the original flags of the Slatec
        routines.
*/

#ifndef BESSEL_LIBRARY_SLATEC_FLAGS_IMPL_H
#define BESSEL_LIBRARY_SLATEC_FLAGS_IMPL_H

#include <stdio.h>

/* Flags for slatec_zbesj_impl_ function */
static inline void slatec_flags_zbesj_impl_(const int ierr, const int nz)
{
    switch(ierr)
    {
        case 0:
            break;
        case 1:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesj_impl_() ->"
                            " Input error: No computation.\n");
            break;
        case 2:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesj_impl_() ->"
                            " Overflow: No computation, imag(z) too large on"
                            " nonscaled calculation.\n");
            break;
        case 3:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesj_impl_() ->"
                            " The abs(z) or abs(nu)+n-1 large: Computation"
                            " done but losses of signifcance by argument"
                            " reduction produce less than half of machine"
                            " accuracy.\n");
            break;
        case 4:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesj_impl_() ->"
                            " The abs(z) or abs(nu)+n-1 too large: No"
                            " computation because of complete losses of"
                            " significance by argument reduction.\n");
            break;
        case 5:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesj_impl_() ->"
                            " Error: No computation, algorithm termination"
                            " condition not met.\n");
            break;
    }
    if (nz != 0) fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                                  " slatec_zbesj_impl_() ->"
                                  " Number of components set to zero due to"
                                  " underflow: %d.\n", nz);
}

/* Flags for slatec_zbesy_impl_ function */
static inline void slatec_flags_zbesy_impl_(const int ierr, const int nz)
{
    switch(ierr)
    {
        case 0:
            break;
        case 1:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesy_impl_() ->"
                            " Input error: No computation.\n");
            break;
        case 2:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesy_impl_() ->"
                            " Overflow: No computation, abs(nu) too large or"
                            " abs(z) is too small or both.\n");
            break;
        case 3:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesy_impl_() ->"
                            " The abs(z) or abs(nu)+n-1 large: Computation"
                            " done but losses of signifcance by argument"
                            " reduction produce less than half of machine"
                            " accuracy.\n");
            break;
        case 4:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesy_impl_() ->"
                            " The abs(z) or abs(nu)+n-1 too large: No"
                            " computation because of complete losses of"
                            " significance by argument reduction.\n");
            break;
        case 5:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesy_impl_() ->"
                            " No computation, algorithm termination condition"
                            " not met.\n");
            break;
    }
    if (nz != 0) fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesy_impl_() ->"
                            " Number of components set to zero due to"
                            " underflow: %d.\n", nz);
}

/* Flags for slatec_zbesh_impl_ function */
static inline void slatec_flags_zbesh_impl_(const int ierr, const int nz)
{
    switch(ierr)
    {
        case 0:
            break;
        case 1:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesh_impl_() ->"
                            " Input error: No computation.\n");
            break;
        case 2:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesh_impl_() ->"
                            " Overflow: No computation, abs(nu) too large or"
                            " abs(z) too small or both.\n");
            break;
        case 3:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesh_impl_() ->"
                            " The abs(z) or abs(nu)+n-1 large: Computation"
                            " done but losses of signifcance by argument"
                            " reduction produce less than half of machine"
                            " accuracy.\n");
            break;
        case 4:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesh_impl_() ->"
                            " The abs(z) or abs(nu)+n-1 too large: No"
                            " computation because of complete losses of"
                            " significance by argument reduction.\n");
            break;
        case 5:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesh_impl_() ->"
                            " No computation, algorithm termination condition"
                            " not met.\n");
            break;
    }
    if (nz != 0) fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                                  " slatec_zbesh_impl_() ->"
                                  " Number of components set to zero due to"
                                  " underflow: %d.\n", nz);
}

/* Flags for slatec_zbesi_impl_ function */
static inline void slatec_flags_zbesi_impl_(const int ierr, const int nz)
{
    switch(ierr)
    {
        case 0:
            break;
        case 1:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesi_impl_() ->"
                            " Input error: No computation.\n");
            break;
        case 2:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesi_impl_() ->"
                            " Overflow: No computation, real(z) too large on"
                            " nonscaled calculation.\n");
            break;
        case 3:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesi_impl_() ->"
                            " The abs(z) or abs(nu)+n-1 large: Computation"
                            " done but losses of signifcance by argument"
                            " reduction produce less than half of machine"
                            " accuracy.\n");
            break;
        case 4:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesi_impl_() ->"
                            " The abs(z) or abs(nu)+n-1 too large: No"
                            " computation because of complete losses of"
                            " significance by argument reduction.\n");
            break;
        case 5:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesi_impl_() ->"
                            " No computation, algorithm termination condition"
                            " not met.\n");
            break;
    }
    if (nz != 0) fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                                  " slatec_zbesi_impl_() ->"
                                  " Number of components set to zero due to"
                                  " underflow: %d.\n", nz);
}

/* Flags for slatec_zbesk_impl_ function */
static inline void slatec_flags_zbesk_impl_(const int ierr, const int nz)
{
    switch(ierr)
    {
        case 0:
            break;
        case 1:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesk_impl_() ->"
                            " Input error: No computation.\n");
            break;
        case 2:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesk_impl_() ->"
                            " Overflow: No computation, abs(nu) is too large"
                            " or abs(z) is too small or both.\n");
            break;
        case 3:
             fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesk_impl_() ->"
                            " The abs(z) or abs(nu)+n-1 large: Computation"
                            " done but losses of signifcance by argument"
                            " reduction produce less than half of machine"
                            " accuracy.\n");
            break;
        case 4:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesk_impl_() ->"
                            " The abs(z) or abs(nu)+n-1 too large: No"
                            " computation because of complete losses of"
                            " significance by argument reduction.\n");
            break;
        case 5:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbesk_impl_() ->"
                            " No computation, algorithm termination condition"
                            " not met.\n");
            break;
    }
    if (nz != 0) fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                                  " slatec_zbesk_impl_() ->"
                                  " Number of components set to zero due to"
                                  " underflow: %d.\n", nz);
}

/* Flags for slatec_zairy_impl_ function */
static inline void slatec_flags_zairy_impl_(const int ierr, const int nz)
{
    switch(ierr)
    {
        case 0:
            break;
        case 1:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zairy_impl_() ->"
                            " Input error: No computation.\n");
            break;
        case 2:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zairy_impl_() ->"
                            " Overflow: No computation,"
                            " real((2/3)*z*sqrt(z)) is too large on"
                            " nonscaled calculation.\n");
            break;
        case 3:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zairy_impl_() ->"
                            " The abs(z) large: Computation completed, losses"
                            " of signifcance by argument reduction produce"
                            " less than half of machine accuracy.\n");
            break;
        case 4:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zairy_impl_() ->"
                            " The abs(z) too large: No computation, complete"
                            " loss of accuracy by argument reduction.\n");
            break;
        case 5:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zairy_impl_() ->"
                            " No computation, algorithm termination condition"
                            " not met.\n");
            break;
    }
    if (nz != 0) fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                                  " slatec_zairy_impl_() ->"
                                  " Number of components set to zero due to"
                                  " underflow: %d.\n", nz);
}

/* Flags for slatec_zbiry_impl_ function */
static inline void slatec_flags_zbiry_impl_(const int ierr)
{
    switch(ierr)
    {
        case 0:
            break;
        case 1:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbiry_impl_() ->"
                            " Input error: No computation.\n");
            break;
        case 2:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbiry_impl_() ->"
                            " Overflow: No computation, real(z) is too large"
                            " on nonscaled calculation.\n");
            break;
        case 3:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbiry_impl_() ->"
                            " The abs(z) large: Computation completed, losses"
                            " of signifcance by argument reduction produce"
                            " less than half of machine accuracy.\n");
            break;
        case 4:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbiry_impl_() ->"
                            " The abs(z) too large: No computation, complete"
                            " loss of accuracy by argument reduction.\n");
            break;
        case 5:
            fprintf(stderr, "[BESSEL-LIBRARY WARNING]"
                            " slatec_zbiry_impl_() ->"
                            " No computation, algorithm termination condition"
                            " not met.\n");
            break;
    }
}

#endif /* BESSEL_LIBRARY_SLATEC_FLAGS_IMPL_H */