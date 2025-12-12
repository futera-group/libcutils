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

#include <cmn/message.h>
#include <cmn/string.h>
#include "mol/distance.h"
#include "mol/icoord.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert xyz struct to zmt molecular data struct
 
   x - pointer to xyz molecular data struct
   t - termplate of z-matrix */
struct zmt_mol *mol_xyz_zmt(struct xyz_mol *x, struct zmt_mol *t) {
  unsigned i,j,id;
  double r,r_min;
  struct zmt_mol *z = NULL;
  z = mol_zmt_new();
  z->n_atoms = x->n_atoms;
  z->atom = mol_zmt_atom_new(z->n_atoms);
  for (i=0; i<z->n_atoms; i++) {
    if (t && (t->atom[i].num!=x->atom[i].num))
      msg_error("atom order in template is inconsistent with input file",1);
    z->atom[i].name = str_copy_new(x->atom[i].name);
    z->atom[i].num = x->atom[i].num;
    /* bonds */
    if (i>0) {
      if (t)
        z->atom[i].bond_id = t->atom[i].bond_id;
      else {
        id = 0;
        r_min = 1.0E+90;
        for (j=0; j<i; j++) {
          r = dist_r(x->atom[i].coord,x->atom[j].coord,3);
          if (r<r_min) {
            r_min = r;
            id = j;
            }
          }
        z->atom[i].bond_id = id;
        }
      z->atom[i].bond_val = icoord_bond(x->atom[i].coord,
        x->atom[z->atom[i].bond_id].coord);
      }
    /* angles */
    if (i>1) {
      if (t)
        z->atom[i].angle_id = t->atom[i].angle_id;
      else {
        if (z->atom[i].bond_id>0 && 
            z->atom[i].bond_id!=z->atom[z->atom[i].bond_id].bond_id)
          z->atom[i].angle_id = (z->atom[z->atom[i].bond_id].bond_id);
        else {
          id = 0;
          r_min = 1.0E+90;
          for (j=0; j<i; j++) {
            if (j!=z->atom[i].bond_id) {
              r = dist_r(x->atom[z->atom[i].bond_id].coord,x->atom[j].coord,3);
              if (r<r_min) {
                r_min = r;
                id = j;
                }
              }
            }
          z->atom[i].angle_id = id;
          }
        }
      z->atom[i].angle_val = icoord_angle(x->atom[i].coord,
        x->atom[z->atom[i].bond_id].coord,
        x->atom[z->atom[i].angle_id].coord);
      }
    /* dihedral angles */
    if (i>2) {
      if (t)
        z->atom[i].dihed_id = t->atom[i].dihed_id;
      else {
        if (z->atom[i].angle_id>0 &&
            z->atom[i].bond_id!=z->atom[z->atom[i].angle_id].bond_id &&
            z->atom[i].angle_id!=z->atom[z->atom[i].angle_id].bond_id)
          z->atom[i].dihed_id = (z->atom[z->atom[i].angle_id].bond_id);
        else {
          id = 0;
          r_min = 1.0E+90;
          for (j=0; j<i; j++) {
            if (j!=z->atom[i].bond_id && j!=z->atom[i].angle_id) {
              r = dist_r(x->atom[z->atom[i].angle_id].coord,x->atom[j].coord,3);
              if (r<r_min) {
                r_min = r;
                id = j;
                }
              }
            }
          z->atom[i].dihed_id = id;
          }
        }
      z->atom[i].dihed_val = icoord_dihed(x->atom[i].coord,
        x->atom[z->atom[i].bond_id].coord,
        x->atom[z->atom[i].angle_id].coord,
        x->atom[z->atom[i].dihed_id].coord);
      }
    }
  return(z);
  }

/* -------------------------------------------------------------------------- */
