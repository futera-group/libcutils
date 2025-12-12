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

/* Add molecular topologies to the system
 
   d - gromacs data structure 
   m - array of molecular data structures
   n - number of molecules */
void gmx_dat_mol_add(struct gmx_dat *d, struct gmx_mol *m, unsigned n) {
  unsigned i,nt;
  struct gmx_mol *t;
  /* new array of molecules */
  nt = d->n_mols + n;
  t = gmx_mol_new(nt);
  for (i=0; i<d->n_mols; i++)
    gmx_mol_copy(t+i,d->mol+i);
  for (i=0; i<n; i++)
    gmx_mol_copy(t+d->n_mols+i,m+i);
  /* replace the array */
  d->mol = gmx_mol_free(d->mol,d->n_mols);
  d->n_mols = nt;
  d->mol = t;
  /* update pointers in fragments */
  for (i=0; i<d->n_frags; i++)
    d->frag[i].mol = gmx_dat_mol_get(d,d->frag[i].name);
  }

/* Return pointer to molecular data structure with specified name
 
   d    - gromacs data structure 
   name - name of the molecule */
struct gmx_mol* gmx_dat_mol_get(struct gmx_dat *d, char *name) {
  unsigned i;
  for (i=0; i<d->n_mols; i++)
    if (str_compare(d->mol[i].name,name))
      return(d->mol+i);
  msg_error_f("molecule '%s' not defined in system topology",1,name);
  return(0);
  }

/* -------------------------------------------------------------------------- */
