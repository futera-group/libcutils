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

/* convert xyz molecular data struct to gen struct

   X - pointer to XYZ data struct */
struct gen_mol* mol_xyz_gen(struct xyz_mol *x) {
  unsigned i,j;
  struct gen_mol *m;
  m = mol_gen_new();
  m->title = str_copy_new(x->title);
  m->n_atoms = x->n_atoms;
  m->atom = mol_gen_atom_new(m->n_atoms);
  for (i=0; i<x->n_atoms; i++) {
    m->atom[i].num = x->atom[i].num;
    for (j=0; j<3; j++)
      m->atom[i].coord[j] = x->atom[i].coord[j];
    }
  return(m);
  }

/* -------------------------------------------------------------------------- */
