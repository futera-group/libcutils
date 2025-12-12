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
#include "mol/cell.h"

/* -------------------------------------------------------------------------- */

/* set cell parameters by given side lengths

   c        - the cell data structure
   a1,a2,a3 - the cell side lengths */
void cell_set_side(struct cell *c, double a1, double a2, double a3) {
  unsigned i,j;
  vec_fset(c->origin,0.0,3);
  vec_fset(c->angle,90.0,3);
  c->side[0] = a1;
  c->side[1] = a2;
  c->side[2] = a3;
  for (i=0; i<3; i++) {
    c->vector[i][i] = c->side[i];
    for (j=i+1; j<3; j++)
      c->vector[i][j] = c->vector[j][i] = 0.0;
    }
  cell_type_set(c);
  }

/* set cell parameters by given array of side lengths

   c - the cell data structure
   a - the side-length array */
void cell_set_side_v(struct cell *c, double *a) {
  cell_set_side(c,a[0],a[1],a[2]);
  }

/* -------------------------------------------------------------------------- */

/* set cell parameters by given side lengths and angles

   c        - the cell data structure
   s1,s2,s3 - the side lengths
   a1,a2,a3 - the angles */
void cell_set_side_angle(struct cell *c, double s1, double s2, double s3,
  double a1, double a2, double a3) {
  double x[3],t[3];
  vec_fset(c->origin,0.0,3);
  c->side[0] = s1;
  c->side[1] = s2;
  c->side[2] = s3;
  c->angle[0] = a1;
  c->angle[1] = a2;
  c->angle[2] = a3;
  cell_tilt_v(c,x,t);
  c->vector[0][0] = x[0];
  c->vector[0][1] = 0.0;
  c->vector[0][2] = 0.0;
  c->vector[1][0] = t[0];
  c->vector[1][1] = x[1];
  c->vector[1][2] = 0.0;
  c->vector[2][0] = t[1];
  c->vector[2][1] = t[2];
  c->vector[2][2] = x[2];
  cell_type_set(c);
  }

/* set cell parameters by given arrays of side lengths and angles

   c - the cell data structure
   s - the side-length array
   a - the angle array */
void cell_set_side_angle_v(struct cell *c, double *s, double *a) {
  cell_set_side_angle(c,s[0],s[1],s[2],a[0],a[1],a[2]);
  }

/* -------------------------------------------------------------------------- */

/* set cell parameters by given side lengths and tilting factors

   c        - the cell data structure
   a1,a2,a3 - the side lengths
   xy,xz,yz - the tilting factors */
void cell_set_side_tilt(struct cell *c, double a1, double a2, double a3,
  double xy, double xz, double yz) {
  double v1[3],v2[3],v3[3];
  /* construct cell vectors */
  v1[0] = a1;
  v1[1] = 0.0;
  v1[2] = 0.0;
  v2[0] = xy;
  v2[1] = a2;
  v2[2] = 0.0;
  v3[0] = xz;
  v3[1] = yz;
  v3[2] = a3;
  /* set the cell parameters */
  cell_set_vec(c,v1,v2,v3);
  cell_type_set(c);
  }

/* set cell parameters by given arrays of side lengths and tilting factors

   c - the cell data structure
   a - the side-length array
   t - the tilting-factor array */
void cell_set_side_tilt_v(struct cell *c, double *a, double *t) {
  cell_set_side_tilt(c,a[0],a[1],a[2],t[0],t[1],t[2]);
  }

/* -------------------------------------------------------------------------- */

/* set cell parameters by given cell vectors

   c        - the cell data structure
   v1,v2,v3 - the cell vectors */
void cell_set_vec(struct cell *c, double *v1, double *v2, double *v3) {
  unsigned i;
  vec_fset(c->origin,0.0,3);
  vec_fcopy(c->vector[0],v1,3);
  vec_fcopy(c->vector[1],v2,3);
  vec_fcopy(c->vector[2],v3,3);
  for (i=0; i<3; i++) {
   c->side[i] = vec_fnorm(c->vector[i],3);
   c->angle[i] = 180.0*vec_fangle(c->vector[i%3],c->vector[(i+1)%3],3)/M_PI;
   }
  cell_type_set(c);
  }

/* -------------------------------------------------------------------------- */
