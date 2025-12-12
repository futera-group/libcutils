/******************************************************************************\
 *                                                                            * 
 *  Libcutils - library of C function                                         * 
 *                                                                            *
 *  Version:             3.4                                                  * 
 *  Date:                28/01/2017                                           *
 *                                                                            * 
 *  Author:              Zdenek Futera                                        * 
 *                                                                            * 
 *  Address:             University College London (UCL)                      * 
 *                       Department of Physics and Astronomy                  * 
 *                       Gower Street, London WC1E 6BT                        * 
 *                       United Kingdom                                       * 
 *                                                                            * 
 *  E-Mail:              z.futera@ucl.ac.uk                                   * 
 *                                                                            * 
\******************************************************************************/

#ifndef ZF_LIB_CMN_ARRAY_H
#define ZF_LIB_CMN_ARRAY_H

/* -------------------------------------------------------------------------- */

/* allocate 2D unsigned integer array with variable lengths */
unsigned **array_u2Valloc(unsigned, unsigned*);
/* allocated 3D unsigned integer array with fixed dimensions */
unsigned ***array_u3alloc(unsigned, unsigned, unsigned);
/* allocated 3D double-precision real array with fixed dimensions */
double ***array_f3alloc(unsigned, unsigned, unsigned);

/* free memory of 2D unsigned integer array with variable lengths */
void* array_u2Vfree(unsigned**, unsigned);
/* free memory reserved for the 3D u-int array with fixed dimensions */
void array_u3free(unsigned***, unsigned, unsigned);
/* free memory reserved for the 3D double array with fixed dimensions */
void array_f3free(double***, unsigned, unsigned);

/* create copy of 2D unsigned integer array with variable lengths */
unsigned **array_u2Vcopy_new(unsigned**, unsigned, unsigned*);

/* search given value in 2D unsigned integer array with variable lengths */
int array_u2Vfind(unsigned**, unsigned, unsigned, unsigned*);

/* change number of rows in array with variable lengths */
unsigned **array_u2Vresize_r(unsigned**, unsigned*, unsigned, unsigned);

/* delete specific row in array with variable lengths */
unsigned **array_u2Vdelete_r(unsigned**, unsigned*, unsigned, unsigned);

/* -------------------------------------------------------------------------- */

#endif
