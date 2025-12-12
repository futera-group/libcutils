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
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* allocate array of vibrational-frequency data structs
 
   n - number of degree-of-freedoms */
struct gauss_freq *gauss_freq_new(unsigned n) {
  unsigned i;
  struct gauss_freq *s = NULL;
  /* memory allocation */
  s = (struct gauss_freq*)malloc(n*sizeof(struct gauss_freq));
  if (!s)
    msg_error("cannot allocate memory for gaussian freq data structure",1);
  /* initialization */
  for (i=0; i<n; i++) {  
    s[i].freq = 0.0;
    s[i].displ = NULL;
    s[i].mu_e = NULL;
    s[i].mu_m = NULL;
    s[i].mu_d = NULL;
    s[i].pol = NULL;
    s[i].rmass = 0.0;
    s[i].fconst = 0.0;
    s[i].ir = 0.0;
    s[i].raman = 0.0;
    s[i].depol_p = 0.0;
    s[i].depol_u = 0.0;
    s[i].strength_dip = 0.0;
    s[i].strength_rot = 0.0;
    s[i].em_angle = 0.0;
    s[i].group = 0;
    }
  return(s);
  }

/* free memory allocated for gaussian frequency data struct array

   g  - pointer to gaussian frequency data struct array
   nf - number of degree-of-freedom
   na - number of atoms */
struct gauss_freq* gauss_freq_free(struct gauss_freq *g, 
  unsigned nf, unsigned na) {
  unsigned i;
  if (g) {
    for (i=0; i<nf; i++) {
      g[i].displ = mat_ffree(g[i].displ,na);
      }
    free(g);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
