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
#include <cmn/vector.h>

/* -------------------------------------------------------------------------- */

/* return dihedral angle among four atoms in degrees

   a1,a2,a3,a4 - cartesian coordinates of the four atoms */
double icoord_dihed(double *a1, double *a2, double *a3, double *a4) {
  short sgn;
  unsigned i;
  double n1[3],n2[3],res;
  double v21[3],v23[3],v32[3],v34[3];
  for (i=0; i<3; i++) {
    v21[i] = a1[i]-a2[i];
    v23[i] = a3[i]-a2[i];
    v32[i] = a2[i]-a3[i];
    v34[i] = a4[i]-a3[i];
    }
  vec_fprod_vector(v21,v23,n1);
  vec_fprod_vector(v32,v34,n2);
  sgn = (vec_fprod_mixed(v23,n1,n2)<0 ? -1 : 1);
  res = vec_fangle(n1,n2,3)*sgn;
  return(res*180.0/M_PI);
  }

/* return dihedral angle and direction among four atoms in degrees

   a1,a2,a3,a4 - cartesian coordinates of the four atoms
   d           - directional vector */
double icoord_dihed_d(double *a1, double *a2, double *a3, double *a4, 
  double *d1, double *d2, double *d3, double *d4) {
  short sgn;
  unsigned i;
  double n1[3],n2[3],res;
  double v21[3],v23[3],v32[3],v34[3];
  for (i=0; i<3; i++) {
    v21[i] = a1[i]-a2[i];
    v23[i] = a3[i]-a2[i];
    v32[i] = a2[i]-a3[i];
    v34[i] = a4[i]-a3[i];
    d1[i] = 0.0; /* DEBUG */
    d2[i] = 0.0; /* DEBUG */
    d3[i] = 0.0; /* DEBUG */
    d4[i] = 0.0; /* DEBUG */
    }
  vec_fprod_vector(v21,v23,n1);
  vec_fprod_vector(v32,v34,n2);
  sgn = (vec_fprod_mixed(v23,n1,n2)<0 ? -1 : 1);
  res = vec_fangle(n1,n2,3)*sgn;
  return(res*180.0/M_PI);
  }

/* -------------------------------------------------------------------------- */
