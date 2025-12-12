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
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* allocate and initialize array of k-point data structs
 
   n - number of the k-points */
struct cp2k_kpoint *cp2k_kpoint_new(unsigned n) {
  unsigned i,j;
  struct cp2k_kpoint *d = NULL;
  /* memory allocation */
  d = (struct cp2k_kpoint*)malloc(n*sizeof(struct cp2k_kpoint));
  if (!d)
    msg_error("cannot allocate memory for cp2k k-point data array",1);
  /* initialization */
  for (i=0; i<n; i++) {
    for (j=0; j<3; j++)
      d[i].crd[j] = 0.0;
    d[i].weight = 0.0;
    d[i].state = NULL;
    }
  return(d);
  }

/* free memory allocated for array of k-point data structs
 
   d  - the array of k-point
   nk - number of the k-points
   ne - number of electronic states (alpha/beta spin)
   ns - number of spin components */
struct cp2k_kpoint* cp2k_kpoint_free(struct cp2k_kpoint *d, unsigned nk, 
  unsigned *ne, unsigned ns) {
  unsigned i;
  if (d) {
    for (i=0; i<nk; i++)
      d[i].state = cp2k_state_free(d[i].state,ne,ns);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
