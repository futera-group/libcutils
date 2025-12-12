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

/* set coordinates to PDB data struct from general molecular struct

   p - pointer to PDB data struct
   m - pointer to general molecular data struct */
void mol_pdb_set_coord_gen(struct pdb_mol *p, struct gen_mol *m) {
  unsigned i,j,k,id=0;
  for (i=0; i<p->n_res; i++)
    for (j=0; j<p->res[i].n_atoms; j++) {
      for (k=0; k<3; k++) 
        p->res[i].atom[j].coord[k] = m->atom[id].coord[k];
      id++;
      }
  }

/* -------------------------------------------------------------------------- */
