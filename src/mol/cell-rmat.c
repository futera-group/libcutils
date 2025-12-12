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
#include <cmn/matrix.h>
#include <cmn/quaternion.h>
#include <cmn/vector.h>
#include "mol/cell.h"

/* -------------------------------------------------------------------------- */

/* calculate axes-rotation matrix aligning A1 to X and placing A2 to XY plane
 
   c - the cell data */
void cell_rmat_set(struct cell *c) {
  double a[3],*q1,*q2,*q,t;
  c->t_rot_a = mat_falloc(3,3);
  /* transformation to standard orientation */
  if (fabs(c->vector[0][1])>CELL_CHECK_ACC ||
      fabs(c->vector[0][2])>CELL_CHECK_ACC ||
      fabs(c->vector[1][2])>CELL_CHECK_ACC) {
    /* align A1 to X axis */
    a[0] = 0.0;
    a[1] = c->vector[0][2];
    a[2] = -c->vector[0][1];
    t = 180.0*acos(c->vector[0][0]/c->side[0])/M_PI;
    q1 = quat_alloc();
    quat_from_axis(q1,a,t);
    /* place A2 to XY plane */
    vec_fcopy(a,c->vector[1],3);
    quat_vec_rot(q1,a,1);
    t = 180.0*acos(a[1]/sqrt(a[1]*a[1]+a[2]*a[2]))/M_PI;
    t = (a[2]<0.0 ? t : -t);
    vec_fcopy(a,c->vector[0],3);
    quat_vec_rot(q1,a,1);
    q2 = quat_alloc();
    quat_from_axis(q2,a,t);
    /* rotation matrix */
    q = quat_alloc();
    quat_mult(q,q2,q1);
    quat_to_mat3(q,c->t_rot_a);
    /* clean memory */
    quat_free(q1);
    quat_free(q2);
    quat_free(q);
    }
  /* axes are already in standard orientation */
  else
    mat_funit(c->t_rot_a,3);
  }

/* -------------------------------------------------------------------------- */
