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

#include <stdio.h>
#include <cmn/string.h>
#include <mol/atom.h>
#include <mol/distance.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* return empirical formula of molecular structure

   d - pointer to gaussian data struct */
char *gauss_mol_formula(struct gauss_dat *d) {
  static char form[1024];
  char *sym,tmp[1024] = "\0";
  unsigned i,id[ATOM_MAX_NUM+1];
  for (i=0; i<(ATOM_MAX_NUM+1); i++)
    id[i] = 0;
  for (i=0; i<d->n_atoms; i++)
    id[d->atom[i].num]++;
  for (i=0; i<(ATOM_MAX_NUM+1); i++)
    if (id[i]) {
      sym = atom_name(i);
      sprintf(form,"%s%s%d",tmp,sym,id[i]);
      sprintf(tmp,"%s",form);
      }
  str_trim(form);
  return(form);
  }

/* -------------------------------------------------------------------------- */

/* calculate nuclear repulsion energy

   d - pointer to gaussian data struct */
double gauss_mol_nucrep(struct gauss_dat *d) {
  double r,en=0.0;
  unsigned i,j;
  for (i=0; i<d->n_atoms; i++)
    for (j=0; j<i; j++) {
      r = dist_r(d->atom[i].coord,d->atom[j].coord,3);
      en += ((d->atom[i].charge[GAUSS_CH_NC]*d->atom[j].charge[GAUSS_CH_NC])/r);
      }
  return(en);
  }

/* -------------------------------------------------------------------------- */

/* get energy of highest occupied molecular orbital
 
   g - pointer to gaussian data struct
   m - molecular orbital ID
   s - spin type (output) */
double gauss_mol_homo(struct gauss_dat *g, unsigned *m, short *s) {
  short spin = GAUSS_SPIN_A;
  unsigned mo_id = 0;
  double en = 0.0; 
  /* spin unrestricted data */
  if (g->mo_b && g->mo_a) {
    if (!g->n_alpha_electrons) {
      en = g->mo_a[g->n_beta_electrons-1].energy;
      mo_id = g->n_beta_electrons-1;
      spin = GAUSS_SPIN_B;
      }
    else if (!g->n_beta_electrons) {
      en = g->mo_a[g->n_alpha_electrons-1].energy;
      mo_id = g->n_alpha_electrons-1;
      spin = GAUSS_SPIN_A;
      }
    else if (g->mo_a[g->n_alpha_electrons-1].energy>=
        g->mo_b[g->n_beta_electrons-1].energy) {
      en = g->mo_a[g->n_alpha_electrons-1].energy;
      mo_id = g->n_alpha_electrons-1;
      spin = GAUSS_SPIN_A;
      }
    else {
      en = g->mo_b[g->n_beta_electrons-1].energy;
      mo_id = g->n_beta_electrons-1;
      spin = GAUSS_SPIN_B;
      }
    }
  /* spin restricted data */
  else if (g->mo_a) {
    en = g->mo_a[g->n_alpha_electrons-1].energy;
    mo_id = g->n_alpha_electrons-1;
    spin = GAUSS_SPIN_A;
    }
  /* finish */
  if (s)
    (*s) = spin;
  if (m)
    (*m) = mo_id;
  return(en);
  }

/* get energy of lowest unoccupied molecular orbital
 
   g - pointer to gaussian data struct
   m - molecular orbital ID
   s - spin type (output) */
double gauss_mol_lumo(struct gauss_dat *g, unsigned *m, short *s) {
  short spin = GAUSS_SPIN_A;
  unsigned mo_id = 0;
  double en = 0.0; 
  /* spin unrestricted data */
  if (g->mo_b && g->mo_a) {
    if (g->mo_a[g->n_alpha_electrons].energy<
        g->mo_b[g->n_beta_electrons].energy) {
      en = g->mo_a[g->n_alpha_electrons].energy;
      mo_id = g->n_alpha_electrons;
      spin = GAUSS_SPIN_A;
      }
    else {
      en = g->mo_b[g->n_beta_electrons].energy;
      mo_id = g->n_beta_electrons;
      spin = GAUSS_SPIN_B;
      }
    }
  /* spin restricted data */
  else if (g->mo_a) {
    en = g->mo_a[g->n_alpha_electrons].energy;
    mo_id = g->n_alpha_electrons;
    spin = GAUSS_SPIN_A;
    }
  /* finish */
  if (s)
    (*s) = spin;
  if (m)
    (*m) = mo_id;
  return(en);
  }

/* -------------------------------------------------------------------------- */
