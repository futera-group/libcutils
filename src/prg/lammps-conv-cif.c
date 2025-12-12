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

#include <stdlib.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include <mol/molec.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* convert LAMMPS data to CIF crystallographic format
 
   d - LAMMPS structure and force-field data */
struct cif_mol *lammps_to_cif(struct lammps_dat *d) {
  unsigned i,j;
  double box[3];
  struct cif_mol *x;
  /* CIF data structure */
  x = mol_cif_new();
  x->n_atoms = d->n_atoms; 
  x->atom = mol_cif_atom_new(x->n_atoms);
  /* cell size */
  cell_set_vec(x->cell,d->box->vector[0],d->box->vector[1],d->box->vector[2]);
  x->cell->space_group_name = str_copy_new("P 1");
  x->cell->space_group_id = 1;
  lammps_box_get_side_v(d->box,box);
  /* atomic data */
  for (i=0; i<x->n_atoms; i++) {
    if (d->n_atom_types && d->type) {
      x->atom[i].num = d->type[d->atom[i].type].num;
      x->atom[i].name = str_copy_new(atom_name(x->atom[i].num));
      }
    else {
      x->atom[i].num = 0;
      x->atom[i].name = str_copy_new("X");
      }
    for (j=0; j<3; j++)
      x->atom[i].coord[j] = (d->atom[i].crd[j] - d->box->min[j])/box[j];
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */
