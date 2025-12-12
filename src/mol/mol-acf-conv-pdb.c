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

/* convert acf molecular data struct to pdb struct
 
   c - pointer to acf molecular data struct */
struct pdb_mol *mol_acf_pdb(struct acf_mol *c) {
  char sym[5];
  unsigned i,j;
  struct pdb_mol *p;
  /* pdb structure */
  p = mol_pdb_new();
  p->title = str_copy_new(c->form);
  p->n_res = 1;
  p->res = mol_pdb_res_new(p->n_res);
  sprintf(sym,"%3s",c->name);
  p->res[0].name = str_copy_new(sym);
  p->res[0].ter = 1;
  p->res[0].id = 1;
  p->res[0].n_atoms = c->n_atoms;
  p->res[0].atom = mol_pdb_atom_new(p->res[0].n_atoms);
  for (i=0; i<c->n_atoms; i++) {
    if (c->atom[i].name && (c->atom[i].name[0]==' ' || 
       (c->atom[i].name[0]>='0' && c->atom[i].name[0]<='9'))) {
      sprintf(sym,"%-4s",c->atom[i].name);
      }
    else {
      sprintf(sym," %-3s",c->atom[i].name);
      }
    p->res[0].atom[i].name = str_copy_new(sym);
    p->res[0].atom[i].charge = c->atom[i].charge;
    p->res[0].atom[i].id = i+1;
    for (j=0; j<3; j++)
      p->res[0].atom[i].coord[j] = c->atom[i].coord[j];
    }
  return(p);
  }

/* -------------------------------------------------------------------------- */
