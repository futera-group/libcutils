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
#include <mol/atom.h>
#include "qmc/basis.h"

/* -------------------------------------------------------------------------- */

/* define user-defined basis set 
 
   b1 - pointer to basis set struct
   b2 - user-defined basis set */
void basis_def_gen(struct basis *b1, struct basis *b2) {
  unsigned ic1,ic2,is,ip;
  struct basis_center *c1,*c2;
  struct basis_shell *s1,*s2;
  /* loop over all atoms */
  for (ic1=0; ic1<b1->n_centers; ic1++) {
    c1 = b1->center+ic1;
    /* locate basis set for the atom */
    c2 = NULL;
    for (ic2=0; ic2<b2->n_centers; ic2++)
      if (c1->type==b2->center[ic2].type) {
        c2 = b2->center+ic2;
        break;
        }
    if (!c2)
      msg_error_f("basis set for %s atom not found\n",1,atom_name(c1->type));
    /* contracted shells */
    c1->n_shells = c2->n_shells;
    c1->shell = basis_shell_new(c1->n_shells);
    for (is=0; is<c1->n_shells; is++) {
      s1 = c1->shell+is;
      s2 = c2->shell+is;
      /* primitive shells */
      s1->type = s2->type;
      s1->n_prim = s2->n_prim;
      s1->exp = vec_falloc(s1->n_prim);
      s1->cf1 = vec_falloc(s1->n_prim);
      if (s1->type==BASIS_SHELL_SP)
        s1->cf2 = vec_falloc(s1->n_prim);
      for (ip=0; ip<s1->n_prim; ip++) {
        s1->exp[ip] = s2->exp[ip];
        s1->cf1[ip] = s2->cf1[ip];
        if (s1->type==BASIS_SHELL_SP)
          s1->cf2[ip] = s1->cf2[ip];
        }
      }
    }
  }

/* -------------------------------------------------------------------------- */
