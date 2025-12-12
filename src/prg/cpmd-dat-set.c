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

#include <stdlib.h>
#include "prg/cpmd.h"

/* -------------------------------------------------------------------------- */

/* set number of types and type ID for each atom
 
   d - pointer to cpmd data struct
   a - atom indicator array */
void cpmd_dat_set_types(struct cpmd_dat *d, unsigned *a) {
  unsigned i,j,id = 0;
  /* number of types */
  d->n_types = 0;
  for (i=0; i<120; i++)
    d->n_types += a[i];
  /* type array */
  d->type = cpmd_type_new(d->n_types);
  for (i=0; i<120; i++) {
    if (a[i]) {
      d->type[id].num = i;
      d->type[id].n_ao_ids = 0;
      d->type[id].n_orbs = 0;
      d->type[id].ao_id = NULL;
      id++;
      }
    }
  /* set type IDs */
  for (i=0; i<d->n_atoms; i++) {
    for (j=0; j<d->n_types; j++)
      if (d->type[j].num==d->atom[i].t_id)
        break;
    d->atom[i].t_id = j;
    }
  }

/* -------------------------------------------------------------------------- */
