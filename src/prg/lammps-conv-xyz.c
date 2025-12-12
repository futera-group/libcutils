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
#include <cmn/vector.h>
#include <mol/atom.h>
#include <mol/molec.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* convert XYZ molecular data to LAMMPS format

   x - the XYZ data */
struct lammps_dat* lammps_from_xyz(struct xyz_mol *x) {
  unsigned i,j;
  short *a;
  struct lammps_dat *d;
  /* memory allocation */
  d = lammps_dat_new();
  /* header */
  d->header = str_copy_new(x->title && x->title[0] ? x->title :
     "# Converted from XYZ");
  /* atomic types */
  a = vec_sialloc(120);
  vec_siset(a,0,120);
  for (i=0; i<x->n_atoms; i++)
    a[x->atom[i].num]++;
  d->n_atom_types = 0;
  for (i=0; i<120; i++)
    if (a[i])
      d->n_atom_types++;
  d->type = lammps_type_new(d->n_atom_types);
  for (i=0,j=0; i<120; i++)
    if (a[i]) {
      d->type[j].num = i;
      d->type[j].name = str_copy_new(atom_name(i));
      d->type[j].mass = atom_mass(i);
      j++;
      }
  vec_sifree(a);
  /* atomic data */
  d->n_atoms = x->n_atoms;
  d->atom = lammps_atom_new(d->n_atoms);
  for (i=0; i<d->n_atoms; i++) {
    for (j=0; j<d->n_atom_types; j++)
      if (x->atom[i].num == d->type[j].num) {
        d->atom[i].type = j;
        break;
        }
    d->atom[i].mol = 1;
    d->atom[i].crd = vec_fcopy_new(x->atom[i].coord,3);
    }
  /* atom style */
  d->atom_style = LAMMPS_ATOM_CHRG;
  /* box */
  if (d->n_atoms) {
    d->box = lammps_box_new();
    d->box->pbc = vec_sialloc(3);
    vec_siset(d->box->pbc,1,3);
    d->box->tilt = vec_falloc(3);
    vec_fset(d->box->tilt,0.0,3);
    d->box->min = vec_falloc(3);
    d->box->max = vec_falloc(3);
    vec_fcopy(d->box->min,d->atom[0].crd,3);
    vec_fcopy(d->box->max,d->atom[0].crd,3);
    for (i=1; i<d->n_atoms; i++)
      for (j=0; j<3; j++) {
        if (d->box->min[j] > d->atom[i].crd[j])
          d->box->min[j] = d->atom[i].crd[j];
        if (d->box->max[j] < d->atom[i].crd[j])
          d->box->max[j] = d->atom[i].crd[j];
        }
    lammps_box_set_vectors(d->box);
    }
  return(d);
  }

/* convert LAMMPS data to XYZ molecular format
 
   d - LAMMPS structure and force-field data */
struct xyz_mol *lammps_to_xyz(struct lammps_dat *d) {
  unsigned i,j;
  struct xyz_mol *x;
  x = mol_xyz_new();
  x->n_atoms = d->n_atoms; 
  x->atom = mol_xyz_atom_new(x->n_atoms);
  for (i=0; i<x->n_atoms; i++) {
    if (d->n_atom_types && d->type)
      x->atom[i].num = d->type[d->atom[i].type].num;
    for (j=0; j<3; j++)
      x->atom[i].coord[j] = d->atom[i].crd[j];
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */
