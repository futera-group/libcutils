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
#include <stdio.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* calculate simulation-box side lengths
 
   t - the box data
   a - the first side (output)
   b - the second side (output)
   c - the third side (output)
   f - open file stream */
void lammps_box_get_side(struct lammps_box *t,
  double *a, double *b, double *c) {
  double lx,ly,lz;
  /* coordinate ranges */
  lx = t->max[0] - t->min[0];
  ly = t->max[1] - t->min[1];
  lz = t->max[2] - t->min[2];
  /* side lengths */
  (*a) = lx;
  (*b) = sqrt(ly*ly + t->tilt[0]*t->tilt[0]);
  (*c) = sqrt(lz*lz + t->tilt[1]*t->tilt[1] + t->tilt[2]*t->tilt[2]);
  }

/* calculate simulation-box side lengths (vector format)
 
   t - the box data
   v - storage for the side lengths
   f - open file stream */
void lammps_box_get_side_v(struct lammps_box *t, double *v) {
  lammps_box_get_side(t,&(v[0]),&(v[1]),&(v[2]));
  }

/* calculate simulation-box angles
 
   t - the box data
   a - the alpha angle (output)
   b - the beta angle (output)
   c - the third angle (output)
   f - open file stream */
void lammps_box_get_angle(struct lammps_box *t,
  double *a, double *b, double *c) {
  double ly,box[3];
  /* coordinate ranges */
  ly = t->max[1] - t->min[1];
  /* box sides */
  lammps_box_get_side_v(t,box);
  /* angles */
  (*a) = 180.0/M_PI*acos(
          (t->tilt[0]*t->tilt[1] + ly*t->tilt[2])/(box[1]*box[2]));
  (*b) = 180.0/M_PI*acos(
          (t->tilt[1]/box[2]));
  (*c) = 180.0/M_PI*acos(
          (t->tilt[0]/box[1]));         
  }

/* calculate simulation-box angles (vector format)
 
   t - the box data
   v - storage for the angles
   f - open file stream */
void lammps_box_get_angle_v(struct lammps_box *t, double *v) {
  lammps_box_get_angle(t,&(v[0]),&(v[1]),&(v[2]));
  }

/* -------------------------------------------------------------------------- */
