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
#include "mol/cell.h"

/* -------------------------------------------------------------------------- */

/* calculate titlting factors of the cell planes

   c        - the cell data structure
   lx,ly,lz - projected side lengts
   xy,xz,yz - the tilting factors */
void cell_tilt(struct cell *c, double *lx, double *ly, double *lz, 
  double *xy, double *xz, double *yz) {
  double s2;
  (*lx) = c->side[0];
  (*xy) = c->side[1]*cos(M_PI*c->angle[2]/180.0);
  (*xz) = c->side[2]*cos(M_PI*c->angle[1]/180.0);
  s2 = pow(c->side[1],2.0)-pow(*xy,2.0);
  (*ly) = (fabs(s2)<1.0e-8 ? 0.0 : sqrt(s2));
  (*yz) = ((c->side[1]*c->side[2]*cos(M_PI*c->angle[0]/180.0)) - 
          ((*xy)*(*xz))) / (*ly);
  s2 = pow(c->side[2],2.0)-pow(*xz,2.0)-pow(*yz,2.0);
  (*lz) = (fabs(s2)<1.0e-8 ? 0.0 : sqrt(s2));
  }

/* calculate arrays of titlting factors of the cell planes

   c - the cell data structure
   s - projected side lengts
   t - the tilting factors */
void cell_tilt_v(struct cell *c, double *s, double *t) {
  cell_tilt(c,&(s[0]),&(s[1]),&(s[2]),&(t[0]),&(t[1]),&(t[2]));
  }

/* -------------------------------------------------------------------------- */
