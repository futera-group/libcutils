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

/* merge two XYZ molecular data structs
 
   x1,x2 - pointer to xyz molecular data structs */
struct xyz_mol *mol_xyz_merge(struct xyz_mol *x1, struct xyz_mol *x2) {
  unsigned i,j,n=0;
  struct xyz_mol *x3;
  /* merged structure */
  x3 = mol_xyz_new();
  x3->n_atoms = (x1->n_atoms+x2->n_atoms);
  x3->atom = mol_xyz_atom_new(x3->n_atoms);
  x3->title = str_copy_new("XYZ merge");
  /* copy first molecule */
  for (i=0; i<x1->n_atoms; i++) {
    x3->atom[n].name = str_copy_new(x1->atom[i].name);
    x3->atom[n].num = x1->atom[i].num;
    for (j=0; j<3; j++)
      x3->atom[n].coord[j] = x1->atom[i].coord[j];
    n++; 
    }
  /* add second molecule */
  for (i=0; i<x2->n_atoms; i++) {
    x3->atom[n].name = str_copy_new(x2->atom[i].name);
    x3->atom[n].num = x2->atom[i].num;
    for (j=0; j<3; j++)
      x3->atom[n].coord[j] = x2->atom[i].coord[j];
    n++;
    }
  return(x3);
  }

/* -------------------------------------------------------------------------- */
