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

/* allocate memory for DL_POLY tethering parameters
   
   n - number of the tethering potentials */
struct dlpoly_prm_teth *dlpoly_prm_teth_new(unsigned n) {
  unsigned i;
  struct dlpoly_prm_teth *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct dlpoly_prm_teth));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY tethering-potential data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].pot = 0;
    d[i].site = 0;
    d[i].v1 = 0.0;
    d[i].v2 = 0.0;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for DL_POLY tethering-potential data 
 
   d - the tethering-potential data */
struct dlpoly_prm_teth *dlpoly_prm_teth_free(struct dlpoly_prm_teth *d) {
  if (d) {
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
