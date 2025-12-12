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

#include <cmn/string.h>
#include <qmc/basis.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* merge data from two gaussian data struct to one 
 
   g - pointer to gaussian data file
   f - name of output file */
struct gauss_dat* gauss_dat_merge(struct gauss_dat *g1, struct gauss_dat *g2) {
  unsigned i;
  struct gauss_dat *d;
  d = gauss_dat_new();
  /* general info */
  d->job_title = str_copy_new("Merged structure");
  d->job_mthd = str_copy_new(str_compare(g1->job_mthd,g2->job_mthd) ? 
    g1->job_mthd : (g1->mo_b || g1->mo_b ? "UX" : "RX"));
  d->job_basis = str_copy_new(str_compare(g1->job_basis,g2->job_basis) ? 
    g1->job_basis : "Gen");
  d->job_mthd_id = (g1->job_mthd_id==g2->job_mthd_id ?
    g1->job_mthd_id : GAUSS_MTHD_UNK);
  d->job_type_id = (g1->job_type_id==g2->job_type_id ? 
    g1->job_type_id : GAUSS_JOB_UNK);
  d->prog_ver = (g1->prog_ver==g2->prog_ver ? 
    g1->prog_ver : GAUSS_PRG_UNK);
  d->prog_rev = (g1->prog_rev==g2->prog_rev ? 
    g1->prog_rev : GAUSS_PRG_UNK);
  /* electrons */
  d->n_electrons = (g1->n_electrons+g2->n_electrons);
  d->n_alpha_electrons = (g1->n_alpha_electrons+g2->n_alpha_electrons);
  d->n_beta_electrons = (g1->n_beta_electrons+g2->n_beta_electrons);
  d->multiplicity = (d->n_alpha_electrons-d->n_beta_electrons+1);
  d->charge = (g1->charge+g2->charge);
  /* atoms */
  d->n_atoms = (g1->n_atoms+g2->n_atoms);
  d->atom = gauss_atom_new(d->n_atoms);
  for (i=0; i<g1->n_atoms; i++)
    gauss_atom_copy(d->atom+i,g1->atom+i);
  for (i=0; i<g2->n_atoms; i++)
    gauss_atom_copy(d->atom+g1->n_atoms+i,g2->atom+i);
  /* basis set */
  d->bs = basis_merge(g1->bs,g2->bs);
  /* molecular orbitals */
  d->mo_a = gauss_mo_merge(
    g1->mo_a,g1->bs->n_ibfce,g1->bs->n_bfce,
    g2->mo_a,g2->bs->n_ibfce,g2->bs->n_bfce);
  d->mo_b = gauss_mo_merge(
    (g1->mo_b ? g1->mo_b : g1->mo_a),g1->bs->n_ibfce,g1->bs->n_bfce,
    (g2->mo_b ? g2->mo_b : g2->mo_a),g2->bs->n_ibfce,g2->bs->n_bfce);
  return(d);
  }

/* -------------------------------------------------------------------------- */
