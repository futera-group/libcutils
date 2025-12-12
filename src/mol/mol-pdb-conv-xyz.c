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

/* convert pdb molecular data struct to xyz struct
 
   p - pointer to pdb molecular data struct */
struct xyz_mol *mol_pdb_xyz(struct pdb_mol *p) {
  unsigned i,j,k,id = 0;
  struct xyz_mol *x;
  /* xyz structure */
  x = mol_xyz_new();
  x->n_atoms = mol_pdb_atom_n(p);
  x->atom = mol_xyz_atom_new(x->n_atoms);
  x->title = str_copy_new(p->title);
  for (i=0; i<p->n_res; i++)
    for (j=0; j<p->res[i].n_atoms; j++) {
      x->atom[id].num = atom_num_pdb(p->res[i].atom[j].name);
      for (k=0; k<3; k++)
        x->atom[id].coord[k] = p->res[i].atom[j].coord[k];
      id++;
      }
  return(x);
  }

/* -------------------------------------------------------------------------- */
