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
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* calculate cell vectors from cell side lengts and tilting factors 
 
   b - the box data */
void lammps_box_set_vectors(struct lammps_box *b) {
  unsigned i,j;
  double side[3],abc[3],ang[3];
  /* initialization */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      b->vector[i][j] = 0.0;
  /* set vectors */
  if (b->pbc[0] && b->pbc[1] && b->pbc[2]) {
    /* side lengths */
    for (i=0; i<3; i++)
      side[i] = b->max[i] - b->min[i];
    /* vector lengths */
    abc[0] = side[0];
    abc[1] = sqrt((side[1]*side[1])+(b->tilt[0]*b->tilt[0]));
    abc[2] = sqrt((side[2]*side[2])+(b->tilt[1]*b->tilt[1])+
                  (b->tilt[2]*b->tilt[2]));
    /* angle cosines */
    ang[0] = ((b->tilt[0]*b->tilt[1])+(side[1]*b->tilt[2]))/(abc[1]*abc[2]);
    ang[1] = b->tilt[1]/abc[2];
    ang[2] = b->tilt[0]/abc[1];
    /* cell vectors */
    b->vector[0][0] = abc[0];
    b->vector[1][0] = abc[1]*ang[2];
    b->vector[1][1] = sqrt((abc[1]*abc[1])-(b->vector[1][0]*b->vector[1][0]));
    b->vector[2][0] = abc[2]*ang[1];
    b->vector[2][1] = ((abc[1]*abc[2]*ang[0])-
                      (b->vector[1][0]*b->vector[2][0]))/b->vector[1][1];
    b->vector[2][2] = sqrt((abc[2]*abc[2])-(b->vector[2][0]*b->vector[2][0])-
                      (b->vector[2][1]*b->vector[2][1]));
    }
  }

/* set box parameters from given side lengths and angle cosines

   t  - the box data
   a  - the first box side length
   b  - the second box side length
   c  - the third box side length
   al - cosine of the alpha angle
   bt - cosine the beta angle
   gm - cosine the gamma angle */
void lammps_box_set(struct lammps_box *t, double a, double b, double c,
  double al, double bt, double gm) {
  double ll[3];
  unsigned i;
  /* tilting factors */
  ll[0] = a;
  t->tilt[0] = b*gm;
  t->tilt[1] = c*bt;
  ll[1] = sqrt(b*b - t->tilt[0]*t->tilt[0]);
  t->tilt[2] = (b*c*al - t->tilt[0]*t->tilt[1])/ll[1];
  ll[2] = sqrt(c*c - t->tilt[1]*t->tilt[1] - t->tilt[2]*t->tilt[2]);
  /* box boundaries */
  for (i=0; i<3; i++)
    t->max[i] = t->min[i] + ll[i];
  /* vectors */
  lammps_box_set_vectors(t);
  }

/* -------------------------------------------------------------------------- */
