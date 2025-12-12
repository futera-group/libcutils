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
#include <cmn/string.h>
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for DL_POLY three-body parameters
   
   n - number of the three-body potentials */
struct dlpoly_prm_tbp *dlpoly_prm_tbp_new(unsigned n) {
  unsigned i;
  struct dlpoly_prm_tbp *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct dlpoly_prm_tbp));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY three-body potential data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].at1 = NULL;
    d[i].at2 = NULL;
    d[i].at3 = NULL;
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

/* free memory allocated for DL_POLY three-body potential data 
 
   d - the three-body potential data
   n - number of the potentials */
struct dlpoly_prm_tbp *dlpoly_prm_tbp_free(struct dlpoly_prm_tbp *d,
  unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++) {
      d[i].at1 = str_free(d[i].at1);
      d[i].at2 = str_free(d[i].at2);
      d[i].at3 = str_free(d[i].at3);
      }
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
