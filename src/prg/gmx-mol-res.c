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
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */


/* Return residuum ID of specified atom in a molecule

   d  - molecular data structure
   ia - ID of the atom in residuum (optional)
   id - ID of the atom in system */
unsigned gmx_mol_res_get_id(struct gmx_mol *d, unsigned *ia, unsigned id) {
  unsigned i,j,is;
  for (i=0,is=0; i<d->n_resids; i++)
    for (j=0; j<d->res[i].n_atoms; j++,is++) {
      if (is==id) {
        if (ia)
          (*ia) = j;
        return(i);
        }
      }
  msg_error_f("atom ID %d not found in molecular data structure",1,id+1);
  return(0);
  }

/* -------------------------------------------------------------------------- */
