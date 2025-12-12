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

/* allocate memory for DL_POLY bond parameters
   
   n - number of the bond potentials */
struct dlpoly_prm_bond *dlpoly_prm_bond_new(unsigned n) {
  unsigned i;
  struct dlpoly_prm_bond *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct dlpoly_prm_bond));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY bond-potential data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].pot = 0;
    d[i].id1 = 0;
    d[i].id2 = 0;
    d[i].v1 = 0.0;
    d[i].v2 = 0.0;
    d[i].v3 = 0.0;
    d[i].v4 = 0.0;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for DL_POLY bond-potential data 
 
   d - the bond-potential data */
struct dlpoly_prm_bond *dlpoly_prm_bond_free(struct dlpoly_prm_bond *d) {
  if (d) {
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
