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

#include <cmn/matrix.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert cif molecular data struct to acf struct
 
   c - pointer to cif molecular data struct */
struct acf_mol *mol_cif_acf(struct cif_mol *c) {
  char sym[128];
  unsigned i,j,at[120];
  struct acf_mol *m;
  /* initialization */
  vec_uset(at,0,120);
  /* acf structure */
  m = mol_acf_new();
  m->n_atoms = c->n_atoms;
  m->atom = mol_acf_atom_new(m->n_atoms);
  /* convertion */
  m->name = str_copy_new("MOL");
  for (i=0; i<m->n_atoms; i++) {
    mol_pdb_atom_name_id(sym,c->atom[i].num,at[c->atom[i].num]++);
    m->atom[i].name = str_copy_new(sym);
    m->atom[i].type = str_copy_new("x");
    if (c->cell->transf)
      mat_fmult_rvec(m->atom[i].coord,c->cell->t_c2r_c,c->atom[i].coord,3,3);
    else
      vec_fcopy(m->atom[i].coord,c->atom[i].coord,3);
    for (j=0; j<3; j++)
      m->atom[i].coord[j] *= c->cell->side[j];
    m->atom[i].num = c->atom[i].num;
    m->atom[i].charge = 0.0;
    }
  /* chemical formula */
  mol_acf_set_formula(m);
  /* interatomic bonds */
  mol_acf_set_bonds(m);
  /* atomic types */
  mol_acf_type_gaff(m);
  return(m);
  }

/* -------------------------------------------------------------------------- */
