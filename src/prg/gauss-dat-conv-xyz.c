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
#include <mol/molec.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* convert gaussian data to xyz molecular file format
 
   d - pointer to gaussian data struct */
struct xyz_mol* gauss_dat_xyz(struct gauss_dat *d) {
  unsigned i,j;
  struct xyz_mol *x;
  x = mol_xyz_new();
  x->n_atoms = d->n_atoms;
  x->atom = mol_xyz_atom_new(x->n_atoms);
  x->title = str_copy_new(d->job_title);
  for (i=0; i<d->n_atoms; i++) {
    x->atom[i].num = d->atom[i].num;
    for (j=0; j<3; j++)
      x->atom[i].coord[j] = d->atom[i].coord[j]*CONV_B_ANG;
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */
