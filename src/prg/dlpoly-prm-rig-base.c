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

/* allocate memory for DL_POLY rigid parameters
   
   n - number of the rigid units */
struct dlpoly_prm_rig *dlpoly_prm_rig_new(unsigned n) {
  unsigned i;
  struct dlpoly_prm_rig *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct dlpoly_prm_rig));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY rigid-unit data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].n_sites = 0;
    d[i].site = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for DL_POLY rigid-unit data 
 
   d - the rigid-unit data
   n - number of the units */
struct dlpoly_prm_rig *dlpoly_prm_rig_free(struct dlpoly_prm_rig *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++)
      if (d[i].site) 
        free(d[i].site);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
