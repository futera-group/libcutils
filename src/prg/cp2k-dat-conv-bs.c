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

#include <cmn/vector.h>
#include <qmc/basis.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* convert cp2k basis set to generic format
 
   d - cp2k data struct */
struct basis *cp2k_dat_conv_basis(struct cp2k_dat *d) {
  unsigned i,j,k;
  struct basis *b;
  struct cp2k_kind *t;
  b = basis_new();
  b->n_centers = d->n_atoms;
  b->center = basis_center_new(b->n_centers);
  for (i=0; i<b->n_centers; i++) {
    t = d->kind+d->atom[i].kind;
    b->center[i].type = d->atom[i].num;
    for (j=0; j<3; j++)
      b->center[i].coord[j] = d->atom[i].crd[j];
    b->center[i].n_shells = t->n_shells;
    b->center[i].shell = basis_shell_new(b->center[i].n_shells);
    for (j=0; j<b->center[i].n_shells; j++) {
      b->center[i].shell[j].type = BASIS_SHELL_X;
      switch (t->shell[j].type) {
        case 0: b->center[i].shell[j].type = BASIS_SHELL_S;  break;
        case 1: b->center[i].shell[j].type = BASIS_SHELL_P;  break;
        case 2: b->center[i].shell[j].type = BASIS_SHELL_Dc; break;
        case 3: b->center[i].shell[j].type = BASIS_SHELL_Fc; break;
        case 4: b->center[i].shell[j].type = BASIS_SHELL_Gc; break;
        case 5: b->center[i].shell[j].type = BASIS_SHELL_Hc; break;
        case 6: b->center[i].shell[j].type = BASIS_SHELL_Ic; break;
        case 7: b->center[i].shell[j].type = BASIS_SHELL_Jc; break;
        }
      if (t->shell[j].n_cart_fce) {
        b->center[i].shell[j].n_prim = t->shell[j].ao[0].n_prim_fce;
        b->center[i].shell[j].exp = vec_falloc(b->center[i].shell[j].n_prim);
        b->center[i].shell[j].cf1 = vec_falloc(b->center[i].shell[j].n_prim);
        for (k=0; k<b->center[i].shell[j].n_prim; k++) {
          b->center[i].shell[j].exp[k] = t->shell[j].ao[0].exp[k];
          b->center[i].shell[j].cf1[k] = t->shell[j].ao[0].coeff[k];
          }
        }
      }
    }
  basis_def_parm(b);
  return(b);
  }

/* -------------------------------------------------------------------------- */
