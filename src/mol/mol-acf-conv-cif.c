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
#include "mol/atom.h"
#include "mol/cell.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert acf struct to cif molecular data struct
 
   c     - pointer to acf molecular data struct
   shift - shift structure to the bottom of the cell
   box_s - box-side lengths
   box_a - box angles */
struct cif_mol *mol_acf_cif(struct acf_mol *c, short shift,
  double *box_s, double *box_a) {
  double max[3] = {-9.0E+90,-9.0E+90,-9.0E+90};
  double min[3] = { 9.0E+90, 9.0E+90, 9.0E+90};
  double crd[3];
  unsigned i,j;
  struct cif_mol *x;
  /* cif structure */
  x = mol_cif_new();
  x->n_atoms = c->n_atoms;
  x->atom = mol_cif_atom_new(x->n_atoms);
  /* smallest / largest coordinate */
  if (shift || !box_s) {
    for (i=0; i<c->n_atoms; i++) {
      for (j=0; j<3; j++) {
        if (c->atom[i].coord[j]<min[j])
          min[j] = c->atom[i].coord[j];
        if (c->atom[i].coord[j]>max[j])
          max[j] = c->atom[i].coord[j];
        }
      }
    }
  /* unit cell */
  if (box_s && box_a)
    cell_set_side_angle_v(x->cell,box_s,box_a);
  else if (box_s)
    cell_set_side_v(x->cell,box_s);
  else
    cell_set_side_angle(x->cell,
      max[0]-min[0],max[1]-min[1],max[2]-min[2],
      90.0,90.0,90.0);
  x->cell->space_group_name = str_copy_new("P 1");
  x->cell->space_group_id = 1;
  /* copy coordinates */
  for (i=0; i<x->n_atoms; i++) {
    x->atom[i].name = str_copy_new(atom_name(c->atom[i].num));
    x->atom[i].num = c->atom[i].num;
    for (j=0; j<3; j++) {
      crd[j] = c->atom[i].coord[j];
      if (shift)
        crd[j] -= min[j];
      crd[j] /= x->cell->side[j];
      }
    if (x->cell->transf) 
      mat_fmult_rvec(x->atom[i].coord,x->cell->t_rot_a,crd,3,3);
    else {
      for (j=0; j<3; j++)
        x->atom[i].coord[j] = crd[j];
      }
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */
