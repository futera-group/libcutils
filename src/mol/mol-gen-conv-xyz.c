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
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert generic molecular data struct to xyz format

   m - pointer to generic data struct */
struct xyz_mol* mol_gen_xyz(struct gen_mol *m) {
  unsigned i,j;
  struct xyz_mol *x;
  x = mol_xyz_new();
  x->title = str_copy_new(m->title);
  x->n_atoms = m->n_atoms;
  x->atom = mol_xyz_atom_new(x->n_atoms);
  for (i=0; i<m->n_atoms; i++) {
    x->atom[i].num = m->atom[i].num;
    for (j=0; j<3; j++)
      x->atom[i].coord[j] = m->atom[i].coord[j];
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */
