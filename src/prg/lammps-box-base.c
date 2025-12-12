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

#include <stdlib.h>
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/vector.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for LAMMPS simulation box data */
struct lammps_box *lammps_box_new(void) {
  struct lammps_box *d = NULL;
  /* memory allocation */
  d = malloc(sizeof(struct lammps_box));
  if (!d) 
    msg_error("cannot allocate memory for LAMMPS simulation box data",1);
  /* initialization */
  d->pbc = vec_sialloc(3);
  vec_siset(d->pbc,0,3);
  d->min = vec_falloc(3);
  vec_fset(d->min,0.0,3);
  d->max = vec_falloc(3);
  vec_fset(d->max,0.0,3);
  d->tilt = vec_falloc(3);
  vec_fset(d->tilt,0.0,3);
  d->vector = mat_falloc(3,3);
  mat_fset(d->vector,0.0,3,3);
  return(d);
  }

/* free memory allocated for LAMMPS simulation box data 
 
   d - the simulation box data */
struct lammps_box *lammps_box_free(struct lammps_box *d) {
  if (d) {
    d->pbc = vec_sifree(d->pbc);
    d->min = vec_ffree(d->min);
    d->max = vec_ffree(d->max);
    d->tilt = vec_ffree(d->tilt);
    d->vector = mat_ffree(d->vector,3);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
