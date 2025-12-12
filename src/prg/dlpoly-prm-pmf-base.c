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

/* allocate memory for DL_POLY PMF parameters
   
   n - number of the PMF potentials */
struct dlpoly_prm_pmf *dlpoly_prm_pmf_new(unsigned n) {
  unsigned i;
  struct dlpoly_prm_pmf *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct dlpoly_prm_pmf));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY PMF data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].id1 = 0;
    d[i].id2 = 0;
    d[i].w1 = 0.0;
    d[i].w2 = 0.0;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for DL_POLY PMF data 
 
   d - the PMF data */
struct dlpoly_prm_pmf *dlpoly_prm_pmf_free(struct dlpoly_prm_pmf *d) {
  if (d) {
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
