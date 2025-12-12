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

/* convert cif struct to xyz molecular data struct

   c - pointer to cif molecular data struct */
struct xyz_mol *mol_cif_xyz(struct cif_mol *c) {
  unsigned i,j;
  struct xyz_mol *x;
  x = mol_xyz_new();
  x->n_atoms = c->n_atoms;
  x->atom = mol_xyz_atom_new(x->n_atoms);
  x->title = str_copy_new("Converted from XYZ");
  for (i=0; i<c->n_atoms; i++) {
    x->atom[i].num = c->atom[i].num;
    if (c->cell->transf)
      mat_fmult_rvec(x->atom[i].coord,c->cell->t_c2r_c,c->atom[i].coord,3,3);
    else
      vec_fcopy(x->atom[i].coord,c->atom[i].coord,3);
    for (j=0; j<3; j++)
      x->atom[i].coord[j] *= c->cell->side[j];
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */
