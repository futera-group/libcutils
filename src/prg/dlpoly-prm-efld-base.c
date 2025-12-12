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

/* allocate memory for DL_POLY external-field parameters
   
   n - number of the external-field potentials */
struct dlpoly_prm_efld *dlpoly_prm_efld_new(unsigned n) {
  unsigned i;
  struct dlpoly_prm_efld *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct dlpoly_prm_efld));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY external-field data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].pot = 0;
    d[i].v1 = 0.0;
    d[i].v2 = 0.0;
    d[i].v3 = 0.0;
    d[i].v4 = 0.0;
    d[i].v5 = 0.0;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for DL_POLY exteranal-field potential data 
 
   d - the external-field data */
struct dlpoly_prm_efld *dlpoly_prm_efld_free(struct dlpoly_prm_efld *d) {
  if (d) {
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
