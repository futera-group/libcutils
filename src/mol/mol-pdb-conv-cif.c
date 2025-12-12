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

/* convert pdb molecular data struct to cif struct
 
   c     - pointer to cif molecular data struct
   shift - shift structure to the bottom of the cell
   box_s - box-side lengths
   box_a - box angles */
struct cif_mol *mol_pdb_cif(struct pdb_mol *p, short shift,
  double *box_s, double *box_a) {
  double max[3] = {-9.0E+90,-9.0E+90,-9.0E+90};
  double min[3] = { 9.0E+90, 9.0E+90, 9.0E+90};
  double crd[3];
  unsigned i,j,k,id = 0;
  struct cif_mol *c;
  /* cif structure */
  c = mol_cif_new();
  c->n_atoms = mol_pdb_atom_n(p);
  c->atom = mol_cif_atom_new(c->n_atoms);
  /* smallest / largest coordinate */
  if (shift || !box_s) {
    for (i=0; i<p->n_res; i++)
      for (j=0; j<p->res[i].n_atoms; j++)
        for (k=0; k<3; k++) {
          if (p->res[i].atom[j].coord[k]<min[k])
            min[k] = p->res[i].atom[j].coord[k];
          if (p->res[i].atom[j].coord[k]>max[k])
            max[k] = p->res[i].atom[j].coord[k];
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
  for (i=0; i<p->n_res; i++)
    for (j=0; j<p->res[i].n_atoms; j++) {
      c->atom[id].num = atom_num_pdb(p->res[i].atom[j].name);
      c->atom[id].name = str_copy_new(p->res[i].atom[j].name);
      for (k=0; k<3; k++) {
        crd[k] = p->res[i].atom[j].coord[k];
        if (shift)
          crd[k] -= min[k];
        crd[k] /= c->cell->side[k];
        }
      if (c->cell->transf)
        mat_fmult_rvec(c->atom[id].coord,c->cell->t_rot_a,crd,3,3);
      else {
        for (k=0; k<3; k++)
          c->atom[id].coord[k] = crd[k];
        }
      id++;
      }
  return(c);
  }

/* -------------------------------------------------------------------------- */
