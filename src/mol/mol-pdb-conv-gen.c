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

#include "mol/atom.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert pdb molecular data struct to gen struct

   p - pointer to PDB data struct */
struct gen_mol* mol_pdb_gen(struct pdb_mol *p) {
  unsigned i,j,k,id = 0;
  struct gen_mol *m;
  m = mol_gen_new();
  m->n_atoms = mol_pdb_atom_n(p);
  m->atom = mol_gen_atom_new(m->n_atoms);
  for (i=0; i<p->n_res; i++)
    for (j=0; j<p->res[i].n_atoms; j++) {
      m->atom[id].num = atom_num_pdb(p->res[i].atom[j].name);
      for (k=0; k<3; k++)
        m->atom[id].coord[k] = p->res[i].atom[j].coord[k];
      id++;
      }
  return(m);
  }

/* -------------------------------------------------------------------------- */
