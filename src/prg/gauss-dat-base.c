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
#include <qmc/basis.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* allocate new gaussian data struct */
struct gauss_dat *gauss_dat_new(void) {
  struct gauss_dat *g = NULL;
  /* memory allocation */
  g = (struct gauss_dat*)malloc(sizeof(struct gauss_dat));
  if (!g)
    msg_error("cannot allocate new gaussian data struct",1);
  /* initialization */
  g->job_title = NULL;
  g->job_mthd = NULL;
  g->job_basis = NULL;
  g->job_mthd_id = 0;
  g->job_type_id = 0;
  g->prog_ver = 0;
  g->prog_rev = 0;
  g->n_atoms = 0;
  g->n_geom = 0;
  g->n_fdeg = 0;
  g->n_vib_modes = 0;
  g->n_states = 0;
  g->charge = 0;
  g->multiplicity = 0;
  g->n_electrons = 0;
  g->n_alpha_electrons = 0;
  g->n_beta_electrons = 0;
  g->n_bg_charges = 0;
  g->bg_chrg = NULL;
  g->atom = NULL;
  g->bs = NULL;
  g->mo_a = NULL;
  g->mo_b = NULL;
  g->state = NULL;
  g->geom = NULL;
  g->nmr = NULL;
  g->freq = NULL;
  vec_liset(g->iop,0,GAUSS_IOP_MAX);
  return(g);
  }

/* free memory allocated for gaussian data struct
   
   g - pointer to gaussian data struct */
struct gauss_dat* gauss_dat_free(struct gauss_dat *g) {
  if (g) {
    g->bg_chrg = mat_ffree(g->bg_chrg,g->n_bg_charges);
    g->atom = gauss_atom_free(g->atom,g->n_atoms);
    if (g->bs) {
      g->mo_a = gauss_mo_free(g->mo_a,g->bs->n_ibfce);
      g->mo_b = gauss_mo_free(g->mo_b,g->bs->n_ibfce);
      g->bs = basis_free(g->bs);
      }
    g->state = gauss_state_free(g->state,g->n_states);
    g->geom = gauss_geom_free(g->geom,g->n_geom,g->n_atoms);
    g->nmr = gauss_nmr_free(g->nmr,g->n_atoms);
    g->freq = gauss_freq_free(g->freq,g->n_vib_modes,g->n_atoms);
    free(g);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
