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
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* allocate and initialize array of atomic orbitals
 
   n - number of the orbitals */
struct cp2k_ao *cp2k_ao_new(unsigned n) {
  unsigned i,j;
  struct cp2k_ao *d = NULL;
  /* memory allocation */
  d = (struct cp2k_ao*)malloc(n*sizeof(struct cp2k_ao));
  if (!d)
    msg_error("cannot allocate memory for cp2k atomic orbital array",1);
  /* initialization */
  for (i=0; i<n; i++) {
    for (j=0; j<3; j++)
      d[i].qnum[j] = 0;
    d[i].n_prim_fce = 0;
    d[i].exp = NULL;
    d[i].coeff = NULL;
    }
  return(d);
  }

/* free memory allocated for cp2k atomic orbital array
 
   d - the atomic orbital array
   n - number of the orbitals */
struct cp2k_ao* cp2k_ao_free(struct cp2k_ao *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++) {
      d[i].exp = vec_ffree(d[i].exp);
      d[i].coeff = vec_ffree(d[i].coeff);
      }
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
