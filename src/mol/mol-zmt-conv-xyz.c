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

#include <math.h>
#include <cmn/string.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert zmt molecular data struct to xyz struct
 
   z - pointer to xyz molecular data struct */
struct xyz_mol *mol_zmt_xyz(struct zmt_mol *z) {
  double v[3],val;
  unsigned i,j;
  struct xyz_mol *x;
  x = mol_xyz_new();
  x->n_atoms = z->n_atoms;
  x->atom = mol_xyz_atom_new(x->n_atoms);
  x->title = str_copy_new("Converted from ZMT");
  /* convert molecular representation */
  for (i=0; i<z->n_atoms; i++) {
    x->atom[i].name = str_copy_new(z->atom[i].name);
    x->atom[i].num = z->atom[i].num;
    /* initialization */
    for (j=0; j<3; j++)
      x->atom[i].coord[j] = v[j] = 0.0;
    /* r-vector */
    if (i==1)
      v[0] = z->atom[i].bond_val;
    else if (i==2) {
      val = (z->atom[i].angle_val<90.0 ?
        180.0-z->atom[i].angle_val : z->atom[i].angle_val);
      v[0] = z->atom[i].bond_val*cos(M_PI*val/180.0);
      v[1] = z->atom[i].bond_val*sin(M_PI*val/180.0);
      }
    else if (i>2) {
      mol_zmt_atom_vec(v,
        x->atom[z->atom[i].bond_id].coord,
        x->atom[z->atom[i].angle_id].coord,
        x->atom[z->atom[i].dihed_id].coord,
        z->atom[i].bond_val,
        z->atom[i].angle_val,
        z->atom[i].dihed_val);
      }
    /* coordinates */
    if (i>0) {
      for (j=0; j<3; j++)
        x->atom[i].coord[j] = x->atom[z->atom[i].bond_id].coord[j]+v[j];
      }
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */
