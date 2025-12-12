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

#include <math.h>
#include "mol/distance.h"

/* -------------------------------------------------------------------------- */

/* return bond distance between two atoms

   a1,a2 - cartesian coordinates of two atoms */
double icoord_bond(double *a1, double *a2) {
  return(dist_r(a1,a2,3));
  }

/* -------------------------------------------------------------------------- */

/* return bond distance and direction between two atoms

   a1,a2 - cartesian coordinates of two atoms
   d1,d2 - directional vector */
double icoord_bond_d(double *a1, double *a2, double *d1, double *d2) {
  unsigned i;
  double d[3],v = 0.0;
  /* bond length */
  for (i=0; i<3; i++) {
    d[i] = a2[i]-a1[i];
    v += d[i]*d[i];
    }
  v = sqrt(v);
  /* derivatives */
  if (v > 0.0)
    for (i=0; i<3; i++) {
      d1[i] = -d[i]/v;
      d2[i] =  d[i]/v;
      }
  return(v);
  }

/* -------------------------------------------------------------------------- */
