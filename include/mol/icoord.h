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

#ifndef ZF_LIB_MOL_ICOORD_H
#define ZF_LIB_MOL_ICOORD_H

/* -------------------------------------------------------------------------- */ 
/* return bond distance between two atoms */
double icoord_bond(double*, double*);
/* return bond distance and direction between two atoms */
double icoord_bond_d(double*, double*, double*, double*);

/* return angle among three atoms */
double icoord_angle(double*, double*, double*);
/* return angle and directionamong three atoms in degrees */
double icoord_angle_d(double*, double*, double*, double*, double*, double*);

/* return dihedral angle among four atoms */
double icoord_dihed(double*, double*, double*, double*);
/* return dihedral angle and direction among four atoms in degrees */
double icoord_dihed_d(double*, double*, double*, double*,
  double*, double*, double*, double*);

/* -------------------------------------------------------------------------- */ 
#endif
