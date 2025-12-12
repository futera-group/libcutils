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
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Add new atoms to residuum data structure
  
   r - the residuum data structure
   a - the atom array
   n - number of atoms */
void gmx_res_atom_add(struct gmx_res *r, struct gmx_atom *a, unsigned n) {
  unsigned i,nt;
  struct gmx_atom *t;
  /* new array of atoms */
  nt = r->n_atoms + n;
  t = gmx_atom_new(nt);
  for (i=0; i<r->n_atoms; i++)
    gmx_atom_copy(t+i,r->atom+i);
  for (i=0; i<n; i++)
    gmx_atom_copy(t+r->n_atoms+i,a+i);
  /* replace the array */
  r->atom = gmx_atom_free(r->atom,r->n_atoms);
  r->n_atoms = nt;
  r->atom = t;
  }

/* Return ID of specified atom in the residuum
 
   d    - the residuum data structure
   name - name of the atom */
unsigned gmx_res_atom_get_id(struct gmx_res *d, char *name) {
  unsigned i;
  for (i=0; i<d->n_atoms; i++)
    if (str_compare(d->atom[i].name,name))
      return(i);
  msg_error_f("atom \"%s\" not found in residuum #%d (\"%s\")",
    1,name,d->id+1,d->name);
  return(0);
  }

/* -------------------------------------------------------------------------- */
