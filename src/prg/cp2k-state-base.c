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
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* allocate and initialize array of state data structs
 
   n - number of the states */
struct cp2k_state *cp2k_state_new(unsigned n) {
  unsigned i;
  struct cp2k_state *d = NULL;
  /* memory allocation */
  d = (struct cp2k_state*)malloc(n*sizeof(struct cp2k_state));
  if (!d)
    msg_error("cannot allocate memory for cp2k state data array",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].energy = 0.0;
    d[i].occup = 0.0;
    d[i].ao_pj = NULL;
    }
  return(d);
  }

/* free memory allocated for array of state data structs
 
   d - the array of states
   n - number of the states
   s - number of spin components */
struct cp2k_state** cp2k_state_free(struct cp2k_state **d, 
  unsigned *n, unsigned s) {
  unsigned i,j;
  if (d) {
    for (i=0; i<s; i++)
      if (d[i]) {
        for (j=0; j<n[i]; j++)
          d[i][j].ao_pj = vec_ffree(d[i][j].ao_pj);
        }
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
