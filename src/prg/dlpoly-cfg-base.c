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
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for DL_POLY structure data */
struct dlpoly_cfg *dlpoly_cfg_new(void) {
  struct dlpoly_cfg *d = NULL;
  /* memory allocation */
  d = malloc(sizeof(struct dlpoly_cfg));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY structure data",1);
  /* initialization */
  d->header = NULL;
  d->data = 0;
  d->pbc = 0; 
  d->n_atoms = 0;
  d->cell = NULL;
  d->sym = NULL;
  d->crd = NULL;
  d->vel = NULL; 
  d->frc = NULL;
  return(d);
  }

/* free memory allocated for DL_POLY structure data 
 
   d - the structure data */
struct dlpoly_cfg *dlpoly_cfg_free(struct dlpoly_cfg *d) {
  if (d) {
    d->header = str_free(d->header);
    d->cell = mat_ffree(d->cell,3);
    d->sym = vec_sfree(d->sym,d->n_atoms);
    d->crd = mat_ffree(d->crd,d->n_atoms);
    d->vel = mat_ffree(d->vel,d->n_atoms);
    d->frc = mat_ffree(d->frc,d->n_atoms);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
