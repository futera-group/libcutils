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
#include <mol/cell.h>
#include "prg/cpmd.h"

/* -------------------------------------------------------------------------- */

/* allocate and initialize new cpmd data struct */
struct cpmd_dat *cpmd_dat_new(void) {
  struct cpmd_dat *d = NULL;
  /* memory allocation */
  d = (struct cpmd_dat*)malloc(sizeof(struct cpmd_dat));
  if (!d)
    msg_error("cannot allocate memory for cpmd data structure",1);
  /* initialization */
  d->n_atoms = 0;
  d->n_types = 0;
  d->n_states = 0;
  d->n_alpha_states = 0;
  d->n_beta_states = 0;
  d->n_electrons = 0;
  d->n_orbitals = 0;
  d->n_geoms = 0;
  d->n_kpoints = 0;
  d->kp_mesh[0] = 0;
  d->kp_mesh[1] = 0;
  d->kp_mesh[2] = 0;
  d->multiplicity = 1;
  d->lsd_calc = 0;
  d->charge = 0.0;
  d->ovrl = NULL;
  d->cell = cell_new();
  d->type = NULL;
  d->state_a = NULL;
  d->state_b = NULL;
  d->kpoint = NULL;
  d->atom = NULL;
  d->g_final = NULL;
  d->g_opt = NULL;
  return(d);
  }

/* free memory allocated for cpmd data struct 
 
   d - pointer to cpmd data struct */
void cpmd_dat_free(struct cpmd_dat *d) {
  if (d) {
    d->ovrl = mat_ffree(d->ovrl,d->n_orbitals);
    d->cell = cell_free(d->cell);
    cpmd_type_free(d->type,d->n_types);
    cpmd_state_free(d->state_a,d->lsd_calc ? d->n_alpha_states : d->n_states);
    cpmd_state_free(d->state_b,d->n_beta_states);
    cpmd_kpoint_free(d->kpoint,d->n_kpoints);
    cpmd_atom_free(d->atom,d->n_atoms);
    cpmd_geom_free(d->g_opt,d->n_geoms,d->n_atoms);
    cpmd_geom_free(d->g_final,1,d->n_atoms);
    free(d);
    }
  }

/* -------------------------------------------------------------------------- */
