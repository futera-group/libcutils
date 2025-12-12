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
#include <qmc/basis.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* return gaussian internal code of basis set shell type
 
   s - the basis set shell type */
short gauss_bs_shell_id(unsigned s) {
  switch (s) {
    case BASIS_SHELL_S:  return(0);
    case BASIS_SHELL_SP: return(-1);
    case BASIS_SHELL_P:  return(1);
    case BASIS_SHELL_Dp: return(-2);
    case BASIS_SHELL_Dc: return(2);
    case BASIS_SHELL_Fp: return(-3);
    case BASIS_SHELL_Fc: return(3);
    case BASIS_SHELL_Gp: return(-4);
    case BASIS_SHELL_Gc: return(4);
    case BASIS_SHELL_Hp: return(-5);
    case BASIS_SHELL_Hc: return(5);
    case BASIS_SHELL_Ip: return(-6);
    case BASIS_SHELL_Ic: return(6);
    case BASIS_SHELL_Jp: return(-7);
    case BASIS_SHELL_Jc: return(7);
    }
  return(0);
  }

/* return basis set shell type ID corresponding to gaussian internal code
 
   s - the gaussian code */
unsigned gauss_bs_shell_type(short s) {
  switch (s) {
    case  0: return(BASIS_SHELL_S);
    case -1: return(BASIS_SHELL_SP);
    case  1: return(BASIS_SHELL_P);
    case -2: return(BASIS_SHELL_Dp);
    case  2: return(BASIS_SHELL_Dc);
    case -3: return(BASIS_SHELL_Fp);
    case  3: return(BASIS_SHELL_Fc);
    case -4: return(BASIS_SHELL_Gp);
    case  4: return(BASIS_SHELL_Gc);
    case -5: return(BASIS_SHELL_Hp);
    case  5: return(BASIS_SHELL_Hc);
    case -6: return(BASIS_SHELL_Ip);
    case  6: return(BASIS_SHELL_Ic);
    case -7: return(BASIS_SHELL_Jp);
    case  7: return(BASIS_SHELL_Jc);
    }
  return(BASIS_SHELL_X);
  }

/* set basis basis set center coodinates according to atoms
 
   g - pointer to gaussian data struct */
void gauss_bs_set_coord(struct gauss_dat *g) {
  unsigned i,j;
  struct basis *b=g->bs;
  if (b->n_centers!=g->n_atoms)
    msg_error("number of basis set centers differs from number of atoms",1);
  for (i=0; i<b->n_centers; i++) {
    if (b->center[i].type!=g->atom[i].num)
      msg_error("basis set center types inconsistent with atomic types",1);
    for (j=0; j<3; j++)
      b->center[i].coord[j] = g->atom[i].coord[j];
    }
  }

/* set basis specified basis set for given structure data
 
   g - pointer to gaussian data struct
   d - user-define basis set */
void gauss_bs_set(struct gauss_dat *g, struct basis *d) {
  unsigned i,j;
  struct basis *b;
  /* clean basis set struct */
  g->bs = basis_free(g->bs);
  g->bs = basis_new();
  b = g->bs;
  /* set centers */
  b->n_centers = g->n_atoms;
  b->center = basis_center_new(b->n_centers);
  for (i=0; i<b->n_centers; i++) {
    b->center[i].type = g->atom[i].num;
    for (j=0; j<3; j++)
      b->center[i].coord[j] = g->atom[i].coord[j];
    }
  /* set basis functions */
  basis_def(b,d,basis_def_id(g->job_basis));
  }

/* -------------------------------------------------------------------------- */
