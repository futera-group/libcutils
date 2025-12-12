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
#include <cmn/vector.h>
#include "mol/atom.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert xyz struct to pdb molecular data struct
 
   x        - pointer to xyz molecular data struct
   std_name - standard atom name without number */
struct pdb_mol *mol_xyz_pdb(struct xyz_mol *x, short std_name) {
  char sym[256];
  unsigned i,j,at[150];
  struct pdb_mol *p;
  vec_uset(at,0,150);
  /* create molecular struct with 1 residuum */
  p = mol_pdb_new();
  p->title = str_copy_new(x->title);
  p->n_res=1;
  p->res = mol_pdb_res_new(1);
  p->res[0].name = str_copy_new("MOL");
  p->res[0].ter = 1;
  p->res[0].id = 1;
  p->res[0].n_atoms = x->n_atoms;
  p->res[0].atom = mol_pdb_atom_new(x->n_atoms);
  /* copy coordinates and set atom names */
  for (i=0; i<x->n_atoms; i++) {
    if (std_name)
      mol_pdb_atom_name_pos(p,0,i,atom_name(x->atom[i].num));
    else {
      if (x->atom[i].name)
        p->res[0].atom[i].name = str_copy_new(x->atom[i].name);
      else {
        mol_pdb_atom_name_id(sym,x->atom[i].num,at[x->atom[i].num]++);
        p->res[0].atom[i].name = str_copy_new(sym);
        }
      }
    p->res[0].atom[i].id = i+1;
    p->res[0].atom[i].charge = 0.0;
    for (j=0; j<3; j++)
      p->res[0].atom[i].coord[j] = x->atom[i].coord[j];
    }
  return(p);
  }

/* -------------------------------------------------------------------------- */
