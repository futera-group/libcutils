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

/* allocate memory for DL_POLY dihedral angle parameters
   
   n - number of the dihedral-angle potentials */
struct dlpoly_prm_dihed *dlpoly_prm_dihed_new(unsigned n) {
  unsigned i;
  struct dlpoly_prm_dihed *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct dlpoly_prm_dihed));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY dihedral-potential data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].pot = 0;
    d[i].id1 = 0;
    d[i].id2 = 0;
    d[i].id3 = 0;
    d[i].id4 = 0;
    d[i].v1 = 0.0;
    d[i].v2 = 0.0;
    d[i].v3 = 0.0;
    d[i].v4 = 0.0;
    d[i].v5 = 0.0;
    d[i].v6 = 0.0;
    d[i].v7 = 0.0;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for DL_POLY dihedral angle-potential data 
 
   d - the dihedral angle-potential data */
struct dlpoly_prm_dihed *dlpoly_prm_dihed_free(struct dlpoly_prm_dihed *d) {
  if (d) {
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
