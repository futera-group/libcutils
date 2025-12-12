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

/* convert apc molecular data struct to cif struct
 
   p     - pointer to amber prep data struct
   shift - shift structure to the bottom of the cell
   box_s - box-side lengths
   box_a - box angles */
struct cif_mol *mol_apc_cif(struct apc_mol *p, short shift,
  double *box_s, double *box_a) {
  double max[3] = {-9.0E+90,-9.0E+90,-9.0E+90};
  double min[3] = { 9.0E+90, 9.0E+90, 9.0E+90};
  double crd[3];
  unsigned i,j;
  struct cif_mol *c;
  /* cif structure */
  c = mol_cif_new();
  c->n_atoms = p->n_atoms;
  c->atom = mol_cif_atom_new(c->n_atoms);
  /* smallest / largest coordinate */
  if (shift || !box_s) {
    for (i=0; i<p->n_atoms; i++) {
      for (j=0; j<3; j++) {
        if (p->atom[i].coord[j]<min[j])
          min[j] = p->atom[i].coord[j];
        if (p->atom[i].coord[j]>max[j])
          max[j] = p->atom[i].coord[j];
        }
      }
    }
  /* unit cell */
  if (box_s && box_a)
    cell_set_side_angle_v(c->cell,box_s,box_a);
  else if (box_s)
    cell_set_side_v(c->cell,box_s);
  else
    cell_set_side_angle(c->cell,
      max[0]-min[0],max[1]-min[1],max[2]-min[2],
      90.0,90.0,90.0);
  c->cell->space_group_name = str_copy_new("P 1");
  c->cell->space_group_id = 1;
  /* copy coordinates */
  for (i=0; i<c->n_atoms; i++) {
    c->atom[i].name = str_copy_new(atom_name_pdb(p->atom[i].name));
    c->atom[i].num = atom_num(c->atom[i].name);
    for (j=0; j<3; j++) {
      crd[j] = p->atom[i].coord[j];
      if (shift)
        crd[j] -= min[j];
      crd[j] /= c->cell->side[j];
      }
    if (c->cell->transf) 
      mat_fmult_rvec(c->atom[i].coord,c->cell->t_rot_a,crd,3,3);
    else {
      for (j=0; j<3; j++)
        c->atom[i].coord[j] = crd[j];
      }
    }
  return(c);
  }

/* -------------------------------------------------------------------------- */
