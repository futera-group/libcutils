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
#include <cmn/vector.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* allocate vector of NMR shift data structs
 
   n - number of states */
struct gauss_nmr *gauss_nmr_new(unsigned n) {
  unsigned i;
  struct gauss_nmr *s = NULL;
  /* memory allocation */
  s = (struct gauss_nmr*)malloc(n*sizeof(struct gauss_nmr));
  if (!s)
    msg_error("cannot allocate memory for gaussian nmr data structure",1);
  /* initialization */
  for (i=0; i<n; i++) {  
    s[i].isotropic = 0.0;
    s[i].anisotropy = 0.0;
    s[i].tensor = mat_falloc(3,3);
    mat_fset(s[i].tensor,0.0,3,3);
    s[i].eigen = vec_falloc(3);
    vec_fset(s[i].eigen,0.0,3);
    }
  return(s);
  }

/* free memory allocated for gaussian nmr data struct

   g - pointer to gaussian nmr data struct
   n - number of atoms */
struct gauss_nmr* gauss_nmr_free(struct gauss_nmr *g, unsigned n) {
  unsigned i;
  if (g) {
    for (i=0; i<n; i++) {
      g[i].tensor = mat_ffree(g[i].tensor,3);
      g[i].eigen = vec_ffree(g[i].eigen);
      }
    free(g);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
