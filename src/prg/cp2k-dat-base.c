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
#include <mol/cell.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* allocate and initialize new cp2k data struct */
struct cp2k_dat *cp2k_dat_new(void) {
  struct cp2k_dat *d = NULL;
  /* allocate memory */
  d = (struct cp2k_dat*)malloc(sizeof(struct cp2k_dat));
  if (!d)
    msg_error("cannot allocate memory for cp2k data struct",1);
  /* initialization */
  d->version = NULL;
  d->project_name = NULL;
  d->run_type = 0;
  d->print_level = 0;
  d->charge = 0;
  d->multiplicity = 0;
  d->n_atom_kinds = 0;
  d->n_atoms = 0;
  d->n_shell_sets = 0;
  d->n_shells = 0;
  d->n_prim_fce = 0;
  d->n_cart_fce = 0;
  d->n_sphr_fce = 0;
  d->n_bfce = 0;
  d->n_ibfce = 0;
  d->n_spins = 0;
  d->n_kpoints = 0;
  d->n_electrons = NULL;
  d->n_states = NULL;
  d->fermi = 0.0;
  d->ovrl = NULL;
  d->tmat = NULL;
  d->cell = cell_new();
  d->atom = NULL;
  d->kind = NULL;
  d->ecpl = NULL;
  d->kpoint = NULL;
  return(d);
  }

/* free memory allocated for cp2k data struct 
 
   d - pointer to cp2k data struct */
struct cp2k_dat* cp2k_dat_free(struct cp2k_dat *d) {
  if (d) {
    d->n_electrons = vec_ufree(d->n_electrons);
    d->version = str_free(d->version);
    d->project_name = str_free(d->project_name);
    d->ovrl = mat_ffree(d->ovrl,d->n_ibfce);
    d->tmat = mat_ffree(d->tmat,d->n_ibfce);
    d->cell = cell_free(d->cell);
    d->atom = cp2k_atom_free(d->atom);
    d->kind = cp2k_kind_free(d->kind);
    d->ecpl = cp2k_ecpl_free(d->ecpl,d->n_spins);
    d->kpoint = cp2k_kpoint_free(d->kpoint,d->n_kpoints,d->n_states,d->n_spins);
    d->n_states = vec_ufree(d->n_states);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
