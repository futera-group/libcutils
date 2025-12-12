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
#include <cmn/vector.h>
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for DL_POLY tersoff parameters
   
   n - number of the tersoff potentials */
struct dlpoly_prm_ters *dlpoly_prm_ters_new(unsigned n) {
  unsigned i;
  struct dlpoly_prm_ters *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct dlpoly_prm_ters));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY tersoff potential data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].at1 = NULL;
    d[i].at2 = NULL;
    d[i].pot = 0;
    d[i].type = 0;
    d[i].v = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for DL_POLY tersoff potential data 
 
   d - the tersoff potential data
   n - number of the potentials */
struct dlpoly_prm_ters *dlpoly_prm_ters_free(struct dlpoly_prm_ters *d,
  unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++) {
      d[i].at1 = str_free(d[i].at1);
      d[i].at2 = str_free(d[i].at2);
      d[i].v = vec_ffree(d[i].v);
      }
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
