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

#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* set coordinates to XYZ data struct from general molecular struct

   x - pointer to XYZ data struct
   m - pointer to general molecular data struct */
void mol_xyz_set_coord_gen(struct xyz_mol *x, struct gen_mol *m) {
  unsigned i,j;
  for (i=0; i<m->n_atoms; i++)
    for (j=0; j<3; j++)
      x->atom[i].coord[j] = m->atom[i].coord[j];
  }

/* -------------------------------------------------------------------------- */
