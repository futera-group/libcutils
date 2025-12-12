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
#include <cmn/matrix.h>
#include <cmn/message.h>
#include "prg/molden.h"

/* -------------------------------------------------------------------------- */

/* allocate array of vibrational-frequency data structs
 
   n - number of degree-of-freedoms */
struct mld_freq *mld_freq_new(unsigned n) {
  unsigned i;
  struct mld_freq *s = NULL;
  /* memory allocation */
  s = (struct mld_freq*)malloc(n*sizeof(struct mld_freq));
  if (!s)
    msg_error("cannot allocate memory for molden freq data structure",1);
  /* initialization */
  for (i=0; i<n; i++) {  
    s[i].freq = 0.0;
    s[i].displ = NULL;
    s[i].ir = 0.0;
    }
  return(s);
  }

/* free memory allocated for molden frequency data struct array

   m  - pointer to molden frequency data struct array
   nf - number of degree-of-freedom
   na - number of atoms */
struct mld_freq* mld_freq_free(struct mld_freq *m, 
  unsigned nf, unsigned na) {
  unsigned i;
  if (m) {
    for (i=0; i<nf; i++) {
      m[i].displ = mat_ffree(m[i].displ,na);
      }
    free(m);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
