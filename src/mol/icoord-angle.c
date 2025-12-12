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

/* return angle among three atoms in degrees

   a1,a2,a3 - cartesian coordinates of the three atoms */
double icoord_angle(double *a1, double *a2, double *a3) {
  unsigned i;
  double a[3],b[3],res;
  for (i=0; i<3; i++) {
    a[i] = a1[i]-a2[i];
    b[i] = a3[i]-a2[i];
    }
  res = vec_fangle(a,b,3);
  return(res*180.0/M_PI);
  }

/* return angle and directionamong three atoms in degrees

   a1,a2,a3 - cartesian coordinates of the three atoms
   d        - directional vector */
double icoord_angle_d(double *a1, double *a2, double *a3,
  double *d1, double *d2, double *d3) {
  unsigned i;
  double v1[3],v2[3],v11,v22,v12,s12,alpha;
  /* valence angle */
  for (i=0; i<3; i++) {
    v1[i] = a1[i]-a2[i];
    v2[i] = a3[i]-a2[i];
    }
  alpha = vec_fangle(v1,v2,3);
  /* derivatives */
  v11 = vec_fprod_scalar(v1,v1,3);
  v12 = vec_fprod_scalar(v1,v2,3);
  v22 = vec_fprod_scalar(v2,v2,3);
  s12 = sqrt((v11*v22)-(v12*v12));
  for (i=0; i<3; i++) {
    d1[i] = (v12*v1[i]-v11*v2[i])/(v11*s12);
    d2[i] = (v11*v22*(v1[i]+v2[i])-v12*(v22*v1[i]+v11*v2[i]))/(v11*v22*s12);
    d3[i] = (v12*v2[i]-v22*v1[i])/(v22*s12);
    }
  return(alpha*180.0/M_PI);
  }

/* -------------------------------------------------------------------------- */
