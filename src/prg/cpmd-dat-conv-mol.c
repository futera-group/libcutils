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
#include <mol/cell.h>
#include <mol/molec.h>
#include "prg/cpmd.h"

/* -------------------------------------------------------------------------- */

/* convert cpmd structure to generic molecular format
 
   d - pointer to cpmd data struct
   s - structure specification (-1 = input, 0 = final, > 0 = g.opt.) */
struct gen_mol *cpmd_dat_conv_mol(struct cpmd_dat *d, long int s) {
  char text[80];
  unsigned i,j;
  struct gen_mol *m = NULL;
  m = mol_gen_new();
  /* cell */
  m->cell = cell_copy_new(d->cell);
  /* input structure */
  if (s<0) {
    m->n_atoms = d->n_atoms;
    m->atom = mol_gen_atom_new(m->n_atoms);
    m->title = str_copy_new("CPMD Input Structure");
    for (i=0; i<d->n_atoms; i++) {
      for (j=0; j<3; j++)
        m->atom[i].coord[j] = d->atom[i].coord[j];
      m->atom[i].num = d->type[d->atom[i].t_id].num;
      }
    }
  /* final structure */
  else if (s==0) {
    if (d->g_final) {
      m->n_atoms = d->n_atoms;
      m->atom = mol_gen_atom_new(m->n_atoms);
      m->title = str_copy_new("CPMD Final Structure");
      for (i=0; i<d->n_atoms; i++) {
        for (j=0; j<3; j++)
          m->atom[i].coord[j] = d->g_final->coord[i][j];
        m->atom[i].num = d->type[d->atom[i].t_id].num;
        }
      }
    }
  /* structure from geometry optimization */
  else {
    if (d->g_opt && s<=d->n_geoms) {
      m->n_atoms = d->n_atoms;
      m->atom = mol_gen_atom_new(m->n_atoms);
      sprintf(text,"CPMD Structure #%ld (E = %.6f H)",
        s,d->g_opt[s-1].energy);
      m->title = str_copy_new(text);
      for (i=0; i<d->n_atoms; i++) {
        for (j=0; j<3; j++)
          m->atom[i].coord[j] = d->g_opt[s-1].coord[i][j];
        m->atom[i].num = d->type[d->atom[i].t_id].num;
        }
      }
    }
  return(m);
  }

/* -------------------------------------------------------------------------- */
