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

#include <cmn/message.h>
#include "qmc/basis.h"

/* -------------------------------------------------------------------------- */

/* locate basis function according to its ID and return points to it
 
   b  - pointer to basis set data struct
   id - ID of the basis function
   pc - pointer to basis set center (output)
   ps - pointer to basis set shell (output)
   pf - ID of basis function in the shell (output) */
void basis_bfce_pointer(struct basis *b, unsigned id, struct basis_center **pc,
  struct basis_shell **ps, unsigned *pf) {
  unsigned ic,is,ib,bf_id;
  struct basis_center *c;
  struct basis_shell *s;
  bf_id = 0;
  for (ic=0; ic<b->n_centers; ic++) {
    c = b->center+ic;
    for (is=0; is<c->n_shells; is++) {
      s = c->shell+is;
      for (ib=0; ib<s->n_bfce; ib++) {
        if (bf_id==id) {
          (*pc) = c;
          (*ps) = s;
          (*pf) = ib;
          return;
          }
        bf_id++;
        }
      }
    }
  msg_error("required basis function not found",1);
  }

/* return number of primitive basis set functions

   b - pointer to basis set data struct */
unsigned basis_bfce_nprim(struct basis *b) {
  unsigned ic,is,n_pfce=0;
  for (ic=0; ic<b->n_centers; ic++)
    for (is=0; is<b->center[ic].n_shells; is++)
      n_pfce += (b->center[ic].shell[is].n_prim*
        basis_shell_nbfce(b->center[ic].shell[is].type));
  return(n_pfce);
  }

/* -------------------------------------------------------------------------- */
