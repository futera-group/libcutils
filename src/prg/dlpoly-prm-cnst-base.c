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
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for DL_POLY constrain parameters
   
   n - number of the contraints */
struct dlpoly_prm_cnst *dlpoly_prm_cnst_new(unsigned n) {
  unsigned i;
  struct dlpoly_prm_cnst *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct dlpoly_prm_cnst));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY constrain data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].id1 = 0;
    d[i].id2 = 0;
    d[i].bl = 0.0;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for DL_POLY constrain-parameter data 
 
   d - the constrain data */
struct dlpoly_prm_cnst *dlpoly_prm_cnst_free(struct dlpoly_prm_cnst *d) {
  if (d) {
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
