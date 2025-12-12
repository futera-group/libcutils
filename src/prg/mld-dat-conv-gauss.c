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

#include "cmn/matrix.h"
#include "cmn/vector.h"
#include "prg/gauss.h"
#include "prg/molden.h"

/* -------------------------------------------------------------------------- */

/* convert molden data to gaussian format
 
   m - pointer to molden data struct */
struct gauss_dat* mld_dat_gauss(struct mld_dat *m) {
  unsigned i;
  struct gauss_dat *d;
  /* gaussian data structure */
  d = gauss_dat_new();
  /* atomic data */
  d->n_atoms = m->n_atoms;
  d->atom = gauss_atom_new(d->n_atoms);
  for (i=0; i<m->n_atoms; i++) {
    vec_fcopy(d->atom[i].coord,m->atom[i].crd,3);
    d->atom[i].num = m->atom[i].num;
    }
  /* vibrational frequencies */
  d->n_vib_modes = m->n_vib_modes;
  d->freq = gauss_freq_new(d->n_vib_modes);
  for (i=0; i<m->n_vib_modes; i++) {
    d->freq[i].freq = m->freq[i].freq;
    d->freq[i].displ = mat_fcopy_new(m->freq[i].displ,m->n_atoms,3);
    d->freq[i].ir = m->freq[i].ir;
    /* electric dipole transition moment */
    d->freq[i].mu_e = vec_falloc(3);
    mld_freq_dip(m->atom,m->freq[i].displ,m->n_atoms,d->freq[i].mu_e);
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
