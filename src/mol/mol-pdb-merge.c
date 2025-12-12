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

/* merge two PDB molecular data structs
 
   p1,p2 - pointers to two pdb data structs */
struct pdb_mol *mol_pdb_merge(struct pdb_mol *p1, struct pdb_mol *p2) {
  unsigned i,j,k,id=0,id_r=1,id_a=1;
  struct pdb_mol *p;
  /* allocate memory */
  p = mol_pdb_new();
  p->n_res = (p1->n_res+p2->n_res);
  p->res = mol_pdb_res_new(p->n_res);
  /* copy first PDB data struct */
  for (i=0; i<p1->n_res; i++) {
    p->res[id].name = str_copy_new(p1->res[i].name);
    p->res[id].ter = p1->res[i].ter;
    p->res[id].id = id_r++;
    p->res[id].n_atoms = p1->res[i].n_atoms;
    p->res[id].atom = mol_pdb_atom_new(p->res[id].n_atoms);
    for (j=0; j<p1->res[i].n_atoms; j++) {
      p->res[id].atom[j].name = str_copy_new(p1->res[id].atom[j].name);
      p->res[id].atom[j].id = id_a++;
      p->res[id].atom[j].charge = p1->res[i].atom[j].charge;
      for (k=0; k<3; k++)
        p->res[id].atom[j].coord[k] = p1->res[i].atom[j].coord[k];
      }
    id++;
    }
  /* add second PDB data struct */
  for (i=0; i<p2->n_res; i++) {
    p->res[id].name = str_copy_new(p2->res[i].name);
    p->res[id].ter = p2->res[i].ter;
    p->res[id].id = id_r++;
    p->res[id].n_atoms = p2->res[i].n_atoms;
    p->res[id].atom = mol_pdb_atom_new(p->res[id].n_atoms);
    for (j=0; j<p2->res[i].n_atoms; j++) {
      p->res[id].atom[j].name = str_copy_new(p2->res[i].atom[j].name);
      p->res[id].atom[j].id = id_a++;
      p->res[id].atom[j].charge = p2->res[i].atom[j].charge;
      for (k=0; k<3; k++)
        p->res[id].atom[j].coord[k] = p2->res[i].atom[j].coord[k];
      }
    id++;
    }
  return(p);
  }

/* -------------------------------------------------------------------------- */
