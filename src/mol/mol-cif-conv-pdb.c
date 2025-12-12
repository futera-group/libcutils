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

/* convert cif molecular data struct to pdb struct
 
   c - pointer to cif molecular data struct */
struct pdb_mol *mol_cif_pdb(struct cif_mol *c) {
  char sym[256];
  unsigned i,j,at[120];
  struct pdb_mol *p;
  /* initialization */
  vec_uset(at,0,120);
  /* pdb structure */
  p = mol_pdb_new();
  p->n_res = 1;
  p->res = mol_pdb_res_new(1);
  p->res[0].name = str_copy_new("MOL");
  p->res[0].ter = 1;
  p->res[0].id = 1;
  p->res[0].n_atoms = c->n_atoms;
  p->res[0].atom = mol_pdb_atom_new(p->res[0].n_atoms);
  for (i=0; i<c->n_atoms; i++) {
    mol_pdb_atom_name_id(sym,c->atom[i].num,at[c->atom[i].num]++);
    p->res[0].atom[i].name = str_copy_new(sym);
    p->res[0].atom[i].id = i+1;
    if (c->cell->transf) 
      mat_fmult_rvec(p->res[0].atom[i].coord,c->cell->t_c2r_c,c->atom[i].coord,3,3);
    else
      vec_fcopy(p->res[0].atom[i].coord,c->atom[i].coord,3);
    for (j=0; j<3; j++)
      p->res[0].atom[i].coord[j] *= c->cell->side[j];
    p->res[0].atom[i].charge = 0.0;
    }
  return(p);
  }

/* -------------------------------------------------------------------------- */
