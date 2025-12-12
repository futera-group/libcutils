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
#include <cmn/units.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert cub molecular data struct to xyz struct
 
   c - pointer to cub molecular data struct */
struct xyz_mol *mol_cub_xyz(struct cub_mol *c) {
  unsigned i,j;
  struct xyz_mol *x;
  /* cube structure */
  x = mol_xyz_new();
  x->n_atoms = c->n_atoms;
  x->atom = mol_xyz_atom_new(x->n_atoms);
  x->title = str_copy_new(c->title);
  for (i=0; i<c->n_atoms; i++) {
    for (j=0; j<3; j++)
      x->atom[i].coord[j] = c->atom[i].coord[j]*CONV_B_ANG;
    x->atom[i].num = c->atom[i].num;
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */
