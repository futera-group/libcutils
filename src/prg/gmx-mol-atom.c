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

#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Add new atoms to specific residuum in molecular structure
  
   m  - the molecular data structure
   ir - residuum ID 
   a  - the atom array
   n  - number of atoms */
void gmx_mol_atom_add(struct gmx_mol *m, unsigned ir, 
  struct gmx_atom *a, unsigned n) {
  unsigned i,j;
  /* add atoms to residuum */
  gmx_res_atom_add(m->res+ir,a,n);
  /* update atom IDs */
  m->n_atoms = 0;
  for (i=0; i<m->n_resids; i++)
    for (j=0; j<m->res[i].n_atoms; j++)
      m->res[i].atom[j].id = m->n_atoms++;
  }

/* -------------------------------------------------------------------------- */
