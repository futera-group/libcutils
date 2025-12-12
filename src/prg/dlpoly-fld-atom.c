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
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* return number of atoms defined in the force field
 
   d - the force-field data */
unsigned dlpoly_fld_atom_n(struct dlpoly_fld *d) {
  unsigned i,j,n = 0;
  for (i=0; i<d->n_molecules; i++)
    for (j=0; j<d->mol[i].n_atoms; j++)
      n += (d->mol[i].n_molecules*d->mol[i].atom[j].n_atoms);
  return(n);
  }

/* find specific atom in the force-field data 
 
   d - the force-field data
   s - name of the atom
   m - molecular ID (output)
   a - atom ID (output) */
short dlpoly_fld_atom_find(struct dlpoly_fld *d, char *s,
  unsigned *m, unsigned *a) {
  short found = 0;
  unsigned i,j;
  for (i=0; i<d->n_molecules; i++)
    for (j=0; j<d->mol[i].n_atoms; j++)
      if (str_compare(d->mol[i].atom[j].name,s)) {
        found = 1;
        (*m) = i;
        (*a) = j;
        return(found);
        }
  return(found);
  }

/* -------------------------------------------------------------------------- */
