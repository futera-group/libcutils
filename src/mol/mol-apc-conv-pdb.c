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

/* convert apc molecular data struct to pdb struct
 
   p - pointer to amber prep data struct */
struct pdb_mol *mol_apc_pdb(struct apc_mol *p) {
  char sym[5];
  unsigned i,j;
  struct pdb_mol *x;
  /* pdb structure */
  x = mol_pdb_new();
  x->title = str_copy_new(p->title);
  x->n_res = 1;
  x->res = mol_pdb_res_new(x->n_res);
  x->res[0].name = str_copy_new(p->resname);
  x->res[0].ter = 1;
  x->res[0].id = 1;
  x->res[0].n_atoms = p->n_atoms;
  x->res[0].atom = mol_pdb_atom_new(x->res[0].n_atoms);
  for (i=0; i<x->res[0].n_atoms; i++) {
    sprintf(sym,"%-3s",p->atom[i].name);
    x->res[0].atom[i].name = str_copy_new(sym);
    x->res[0].atom[i].id = i+1;
    x->res[0].atom[i].charge = p->atom[i].charge;
    for (j=0; j<3; j++)
      x->res[0].atom[i].coord[j] = p->atom[i].coord[j];
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */
