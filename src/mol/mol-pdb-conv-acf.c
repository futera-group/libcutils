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

/* convert pdb struct to acf molecular data struct
 
   p - pointer to pdb molecular data struct */
struct acf_mol *mol_pdb_acf(struct pdb_mol *p) {
  unsigned i,j,k,id = 0;
  struct acf_mol *c;
  /* acf structure */
  c = mol_acf_new();
  c->n_atoms = mol_pdb_atom_n(p);
  c->atom = mol_acf_atom_new(c->n_atoms);
  /* convert molecular representation */
  c->name = str_copy_new(p->res[0].name);
  for (i=0; i<p->n_res; i++) {
    for (j=0; j<p->res[i].n_atoms; j++) {
      c->atom[id].name = str_copy_new(p->res[i].atom[j].name);
      c->atom[id].type = str_copy_new("x");
      for (k=0; k<3; k++)
         c->atom[id].coord[k] = p->res[i].atom[j].coord[k];
      c->atom[id].num = atom_num_pdb(p->res[i].atom[j].name);
      c->atom[id].charge = p->res[i].atom[j].charge;
      c->charge += (c->atom[id].charge);
      id++;
      }
    }
  /* chemical formula */
  mol_acf_set_formula(c);
  /* interatomic bonds */
  mol_acf_set_bonds(c);
  /* atomic types */
  mol_acf_type_gaff(c);
  return(c);
  }

/* -------------------------------------------------------------------------- */
