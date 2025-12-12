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
#include "mol/atom.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert apc molecular data struct to xyz struct
 
   p - pointer to amber prep data struct */
struct xyz_mol *mol_apc_xyz(struct apc_mol *p) {
  unsigned i,j;
  struct xyz_mol *x;
  /* xyz structure */
  x = mol_xyz_new();
  x->n_atoms = p->n_atoms;
  x->atom = mol_xyz_atom_new(x->n_atoms);
  x->title = str_copy_new(p->title);
  for (i=0; i<x->n_atoms; i++) {
    for (j=0; j<3; j++)
      x->atom[i].coord[j] = p->atom[i].coord[j];
    x->atom[i].num = atom_num_pdb(p->atom[i].name);
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */
