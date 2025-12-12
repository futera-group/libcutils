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

/* convert acf struct to xyz molecular data struct
 
   c - pointer to acf molecular data struct */
struct xyz_mol *mol_acf_xyz(struct acf_mol *c) {
  unsigned i,j;
  struct xyz_mol *x;
  x = mol_xyz_new();
  x->n_atoms = c->n_atoms;
  x->atom = mol_xyz_atom_new(x->n_atoms);
  for (i=0; i<c->n_atoms; i++) {
    for (j=0; j<3; j++)
      x->atom[i].coord[j] = c->atom[i].coord[j];
    x->atom[i].num = c->atom[i].num;
    }
  x->title = str_new(str_length(c->name)+str_length(c->form)+2);
  sprintf(x->title,"%s: %s",c->name,c->form);
  return(x);
  }

/* -------------------------------------------------------------------------- */
