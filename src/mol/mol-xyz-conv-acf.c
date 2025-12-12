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
#include <cmn/vector.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert xyz molecular data struct to acf struct
 
   x - pointer to xyz molecular data struct */
struct acf_mol *mol_xyz_acf(struct xyz_mol *x) {
  char name[256];
  unsigned i,j,at[120];
  struct acf_mol *c;
  c = mol_acf_new();
  c->n_atoms = x->n_atoms;
  c->atom = mol_acf_atom_new(c->n_atoms);
  /* auxiliary array */
  vec_uset(at,0,120);
  /* convert molecular representation */
  c->name = str_copy_new("MOL");
  for (i=0; i<x->n_atoms; i++) {
    if (x->atom[i].name)
      c->atom[i].name = str_copy_new(x->atom[i].name);
    else {
      mol_pdb_atom_name_id(name,x->atom[i].num,at[x->atom[i].num]++);
      c->atom[i].name = str_copy_new(name);
      }
    c->atom[i].type = str_copy_new("x");
    for (j=0; j<3; j++)
      c->atom[i].coord[j] = x->atom[i].coord[j];
    c->atom[i].num = x->atom[i].num;
    c->atom[i].charge = 0.0;
    }
  c->charge = 0.0;
  /* chemical formula */
  mol_acf_set_formula(c);
  /* interatomic bonds */
  mol_acf_set_bonds(c);
  /* atomic types */
  mol_acf_type_gaff(c);
  return(c);
  }

/* -------------------------------------------------------------------------- */
