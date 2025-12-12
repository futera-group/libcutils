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
#include <cmn/message.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for LAMMPS trajectory data */
struct lammps_trj *lammps_trj_new(void) {
  struct lammps_trj *d = NULL;
  /* memory allocation */
  d = malloc(sizeof(struct lammps_trj));
  if (!d) 
    msg_error("cannot allocate memory for LAMMPS trajectory data",1);
  /* initialization */
  d->title[0][0] = '\0';
  d->title[1][0] = '\0';
  d->has_box = 0;
  d->has_4d = 0;
  d->has_chrg = 0;
  d->n_frames = 0;
  d->n_steps = 0;
  d->save_freq = 0;
  d->timestep = 0.0;
  return(d);
  }

/* free memory allocated for LAMMPS trajectory data 
 
   d - the trajectory data */
struct lammps_trj *lammps_trj_free(struct lammps_trj *d) {
  if (d) {
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
