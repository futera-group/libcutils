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
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for new GTO pair array
 
   n - number of the pairs */
struct gto_pair* gto_pair_new(unsigned n) {
  unsigned i,j;
  struct gto_pair *p = NULL;
  /* memory allocation */
  p = (struct gto_pair*)malloc(n*sizeof(struct gto_pair));
  if (!p)
    msg_error("cannot allocate memory for GTO pair array",1);
  /* initialization */
  p->am1 = 0;
  p->am2 = 0;
  p->e1 = 0.0;
  p->e2 = 0.0;
  p->e12 = 0.0;
  p->sigma = 0.0;
  p->pref =0.0;
  for (i=0; i<3; i++) {
    p->aa1[i] = 0;
    p->aa2[i] = 0;
    p->r12[i] = 0.0;
    p->crd[i] = 0.0;
    }
  p->n_pq = 0;
  for (i=0; i<GTO_MAX_I2EXP; i++) {
    for (j=0; j<3; j++)
      p->pq[i].id[j] = 0;
    p->pq[i].coeff = 0.0;
    }
  return(p);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for GTO pair array
 
   p - pointer to pair array */
struct gto_pair* gto_pair_free(struct gto_pair* p) {
  if (p)
    free(p);
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
