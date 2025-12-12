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
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* find atom in specified residuum and returns its ID
  
   r  - the residuum data structure
   a  - name of the atom
   ok - status (optional) */
unsigned mol_pdb_res_atom_find(struct pdb_res *r, char *a, short *ok) {
  unsigned i,id = 0;
  if (ok)
    (*ok) = 0;
  for (i=0; i<r->n_atoms; i++) {
    if (str_compare(r->atom[i].name,a)) {
      if (ok)
        (*ok) = 1;
      id = i;
      break;
      }
    }
  return(id);
  }

/* -------------------------------------------------------------------------- */
