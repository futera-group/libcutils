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

/* allocate vector of state gaussian data struct
 
   n - number of states */
struct gauss_state *gauss_state_new(unsigned n) {
  unsigned i;
  struct gauss_state *s = NULL;
  /* memory allocation */
  s = (struct gauss_state*)malloc(n*sizeof(struct gauss_state));
  if (!s)
    msg_error("cannot allocate memory for gaussian state data structure",1);
  /* intialialization */
  for (i=0; i<n; i++) {
    s[i].energy = 0.0;
    s[i].dipole = NULL;
    s[i].dipole_ee = NULL;
    s[i].lambda = 0.0;
    s[i].strength = 0.0;
    s[i].spin_s2 = 0.0;
    s[i].rot_vel = 0.0;
    s[i].rot_len = 0.0;
    }
  return(s);
  }

/* free memory allocated for gaussian state data struct

   g - pointer to vector of gaussian state data structs 
   n - number of states */
struct gauss_state* gauss_state_free(struct gauss_state *g, unsigned n) {
  unsigned i,j;
  if (g) {
    for (i=0; i<n; i++) {
      g[i].dipole = vec_ffree(g[i].dipole);
      if (g[i].dipole_ee) {
        for (j=0; j<i; j++)
          g[i].dipole_ee[j] = vec_ffree(g[i].dipole_ee[j]);
        free(g[i].dipole_ee);
        }
      }
    free(g);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
