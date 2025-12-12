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

/* allocate array of state data structs

   n - number of states */
struct cpmd_state *cpmd_state_new(unsigned n) {
  unsigned i;
  struct cpmd_state *d = NULL;
  /* memory allocation */
  d = (struct cpmd_state*)malloc(n*sizeof(struct cpmd_state));
  if (!d)
    msg_error("cannot allocate memory for cpmd state data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].energy = 0.0;
    d[i].compln = 0.0;
    d[i].occup = 0.0;
    d[i].ao_pj = NULL;
    }
  return(d);
  }

/* free memory allocated for array of state data structs
 
   d - pointer to the array
   n - number of the structs */
struct cpmd_state *cpmd_state_free(struct cpmd_state *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++)
      d[i].ao_pj = vec_ffree(d[i].ao_pj);
    free(d);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
