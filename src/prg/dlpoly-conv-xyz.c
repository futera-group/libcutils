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

#include <mol/atom.h>
#include <mol/molec.h>
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* convert DL_POLY structure & force-field data to XYZ molecule format
 
   c - the structure data
   f - the force-field data */
struct xyz_mol *dlpoly_to_xyz(struct dlpoly_cfg *c, struct dlpoly_fld *f) {
  unsigned i,j,k,l,m,id;
  struct xyz_mol *x;
  x = mol_xyz_new();
  x->n_atoms = c->n_atoms;
  x->atom = mol_xyz_atom_new(x->n_atoms);
  id = 0;
  for (i=0; i<f->n_molecules; i++)
    for (j=0; j<f->mol[i].n_molecules; j++)
      for (k=0; k<f->mol[i].n_atoms; k++)
        for (l=0; l<f->mol[i].atom[k].n_atoms; l++) {
          x->atom[id].num = atom_num_mass(f->mol[i].atom[k].mass);
          for (m=0; m<3; m++)
            x->atom[id].coord[m] = c->crd[id][m];
          id++;
          }
  return(x);
  }

/* -------------------------------------------------------------------------- */
