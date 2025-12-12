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
#include <cmn/vector.h>
#include "prg/cpmd.h"

/* -------------------------------------------------------------------------- */

/* allocate array of type data structs

   n - number of types */
struct cpmd_type *cpmd_type_new(unsigned n) {
  unsigned i;
  struct cpmd_type *d=NULL;
  /* memory allocation */
  d = (struct cpmd_type*)malloc(n*sizeof(struct cpmd_type));
  if (!d)
    msg_error("cannot allocate memory for cpmd type data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].num = 0;
    d[i].n_ao_ids = 0;
    d[i].n_orbs = 0;
    d[i].charge = 0.0;
    d[i].ao_id = NULL;
    }
  return(d);
  }

/* free memory allocated for array of type data structs
 
   d - pointer to the array
   n - number of the structs */
struct cpmd_type *cpmd_type_free(struct cpmd_type *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++)
      d[i].ao_id = vec_sifree(d[i].ao_id);
    free(d);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
