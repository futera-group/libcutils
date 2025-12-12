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

/* allocate array of k-point data structs

   n - number of k-points */
struct cpmd_kpoint *cpmd_kpoint_new(unsigned n) {
  unsigned i,j;
  struct cpmd_kpoint *d = NULL;
  /* memory allocation */
  d = (struct cpmd_kpoint*)malloc(n*sizeof(struct cpmd_kpoint));
  if (!d)
    msg_error("cannot allocate memory for cpmd kpoint data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    for (j=0; j<3; j++)
      d[i].coord[j] = 0.0;
    d[i].weight = 0.0;
    }
  return(d);
  }

/* free memory allocated for array of k-point data structs
 
   d - pointer to the array
   n - number of the structs */
struct cpmd_kpoint *cpmd_kpoint_free(struct cpmd_kpoint *d, unsigned n) {
  if (d)
    free(d);
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
