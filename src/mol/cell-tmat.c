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
#include <cmn/vector.h>
#include "mol/cell.h"

/* -------------------------------------------------------------------------- */

/* check if the coordinate transformation is needed and set the matrices
 
   c - the cell data */
void cell_tmat_set(struct cell *c) {
  unsigned i,j;
  /* set coordinate-transformation indicator */
  c->transf = 0;
  for (i=0; i<3 && !c->transf; i++) 
    for (j=0; j<3; j++)
      if (i!=j && fabs(c->vector[i][j])>CELL_CHECK_ACC) {
        c->transf = 1;
        break;
        }
  /* create transformation matrices */
  if (c->transf) {
    /* axes: cell -> cartesian */
    c->t_c2r_a = mat_falloc(3,3);
    for (i=0; i<3; i++)
      for (j=0; j<3; j++)
        c->t_c2r_a[i][j] =  c->vector[i][j]/c->side[i];
    /* axes: cartesian -> cell */
    c->t_r2c_a = mat_falloc(3,3);
    mat_finv(c->t_r2c_a,c->t_c2r_a,3);
    /* coord: cell -> cartesian */
    c->t_c2r_c = mat_falloc(3,3);
    for (i=0; i<3; i++)
      for (j=0; j<3; j++)
        c->t_c2r_c[j][i] =  c->vector[i][j]/c->side[i];
    /* coord: cartesian -> cell */
    c->t_r2c_c = mat_falloc(3,3);
    mat_finv(c->t_r2c_c,c->t_c2r_c,3);
    /* axes rotation */
    cell_rmat_set(c);
    }
  }

/* -------------------------------------------------------------------------- */
