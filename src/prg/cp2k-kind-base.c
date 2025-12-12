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

/* allocate and initialize array of atom kinds
 
   n - number of the kinds */
struct cp2k_kind *cp2k_kind_new(unsigned n) {
  unsigned i;
  struct cp2k_kind *d = NULL;
  /* memory allocation */
  d = (struct cp2k_kind*)malloc(n*sizeof(struct cp2k_kind));
  if (!d)
    msg_error("cannot allocate memory for cp2k atom kind array",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].n_atoms = 0;
    d[i].n_shell_sets = 0;
    d[i].n_shells = 0;
    d[i].n_prim_fce = 0;
    d[i].n_cart_fce = 0;
    d[i].n_sphr_fce = 0;
    d[i].norm_type = 0;
    d[i].shell = NULL;
    }
  return(d);
  }

/* free memory allocated for cp2k atom kind array
 
   d - the atom kind array */
struct cp2k_kind* cp2k_kind_free(struct cp2k_kind *d) {
  if (d) {
    d->shell = cp2k_shell_free(d->shell,d->n_shells);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
