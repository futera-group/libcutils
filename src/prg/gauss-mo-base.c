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
#include <cmn/vector.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* allocate vector of state gaussian MO struct
 
   n - number of states */
struct gauss_mo *gauss_mo_new(unsigned n) {
  unsigned i;
  struct gauss_mo *s = NULL;
  /* memory allocation */
  s = (struct gauss_mo*)malloc(n*sizeof(struct gauss_mo));
  if (!s)
    msg_error("cannot allocate memory for gaussion MO data structure",1);
  /* initialization */
  for (i=0; i<n; i++) {
    s[i].energy = 0.0;
    s[i].coeff = NULL;
    }
  return(s);
  }

/* free memory allocated for gaussian MO data struct

   g - pointer to vector of gaussian MO data structs 
   n - number of MOs */
struct gauss_mo* gauss_mo_free(struct gauss_mo *g, unsigned n) {
  unsigned i;
  if (g) {
    for (i=0; i<n; i++)
      g[i].coeff = vec_ffree(g[i].coeff);
    free(g);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
