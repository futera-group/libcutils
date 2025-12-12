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
#include <cmn/message.h>
#include "prg/cpmd.h"

/* -------------------------------------------------------------------------- */

/* allocate array of atom data structs

   n - number of atoms */
struct cpmd_atom *cpmd_atom_new(unsigned n) {
  unsigned i,j;
  struct cpmd_atom *d = NULL;
  /* memory allocation */
  d = (struct cpmd_atom*)malloc(n*sizeof(struct cpmd_atom));
  if (!d)
    msg_error("cannot allocate memory for cpmd atom data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].t_id = 0;
    for (j=0; j<3; j++) {
      d[i].charge[j] = 0.0;
      d[i].coord[j] = 0.0;
      }
    d[i].valence = 0.0;
    }
  return(d);
  }

/* free memory allocated for array of atom data structs
 
   d - pointer to the array
   n - number of the structs */
struct cpmd_atom *cpmd_atom_free(struct cpmd_atom *d, unsigned n) {
  if (d) {
    free(d);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
