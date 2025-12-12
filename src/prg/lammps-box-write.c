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

/* write simulation box parameters to LAMMPS data file 
 
   b - the box data
   f - open file stream */
void lammps_box_write(struct lammps_box *b, FILE *f) {
  /* simulation box */
  if (b->pbc[0])
    fprintf(f,"%12.6f %12.6f xlo xhi\n",b->min[0],b->max[0]);
  if (b->pbc[1])
    fprintf(f,"%12.6f %12.6f ylo yhi\n",b->min[1],b->max[1]);
  if (b->pbc[2])
    fprintf(f,"%12.6f %12.6f zlo zhi\n",b->min[2],b->max[2]);
  if (fabs(b->tilt[0])>1.0e-3 || 
      fabs(b->tilt[1])>1.0e-3 || 
      fabs(b->tilt[2])>1.0e-3)
    fprintf(f,"%12.6f %12.6f %12.6f xy xz yz\n",
      b->tilt[0],b->tilt[1],b->tilt[2]);
  }

/* -------------------------------------------------------------------------- */
