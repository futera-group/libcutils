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
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* allocate and initialize array of orbital shells
 
   n - number of the shells */
struct cp2k_shell *cp2k_shell_new(unsigned n) {
  unsigned i;
  struct cp2k_shell *d = NULL;
  /* memory allocation */
  d = (struct cp2k_shell*)malloc(n*sizeof(struct cp2k_shell));
  if (!d)
    msg_error("cannot allocate memory for cp2k orbital shell array",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].set = 0;
    d[i].type = 0;
    d[i].n_cart_fce = 0;
    d[i].n_sphr_fce = 0;
    d[i].ao = NULL;
    }
  return(d);
  }

/* free memory allocated for cp2k orbital shell array
 
   d - the orbital shell array
   n - number of the shells */
struct cp2k_shell* cp2k_shell_free(struct cp2k_shell *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++)
      d[i].ao = cp2k_ao_free(d[i].ao,d[i].n_cart_fce);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
