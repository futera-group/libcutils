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

#include <stdio.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include "qmc/basis.h"

/* -------------------------------------------------------------------------- */

/* convert basis set type ID to string name

   id - id of the shell */
char* basis_def_name(short id) {
  static char name[10]="\0";
  switch (id) {
    case BASIS_DEF_STO3G: sprintf(name,"STO-3G"); break;
    default:              sprintf(name,"Gen");    break;
    }
  return(name);
  }

/* -------------------------------------------------------------------------- */

/* return ID of basis set type 

   name - name of the shell */
short basis_def_id(char *name) {
  unsigned id = BASIS_DEF_GEN;
  if (str_compare(name,"Gen"))
    id = BASIS_DEF_GEN;
  else if (str_compare(name,"STO3G"))
    id = BASIS_DEF_STO3G;
  else
    msg_error_f("unknown basis set name \"%s\"",1,name);
  return(id);
  }

/* -------------------------------------------------------------------------- */

/* define global parameters of basis set
 
   b - pointer to basis set struct */
void basis_def_parm(struct basis *b) {
  unsigned ic,is,ang_mom,bf_id = 0;
  struct basis_center *c;
  struct basis_shell *s;
  /* initialization */
  b->n_bfce = 0;
  b->n_cont_shells = 0;
  b->n_prim_shells = 0;
  b->n_pure_d = 0;
  b->n_pure_f = 0;
  b->max_ang_mom = 0;
  b->max_bfce_cont = 0;
  /* count functions */
  for (ic=0; ic<b->n_centers; ic++) {
    c = b->center+ic;
    c->bf = bf_id;
    c->n_bfce = 0;
    b->n_cont_shells += c->n_shells;
    for (is=0; is<c->n_shells; is++) {
      s = c->shell+is;
      s->bf = bf_id;
      if (s->type==BASIS_SHELL_Dp)
        b->n_pure_d++;
      if (s->type==BASIS_SHELL_Fp)
        b->n_pure_f++;
      if (s->n_prim>b->max_bfce_cont)
        b->max_bfce_cont = s->n_prim;
      ang_mom = basis_shell_ang_mom(s->type);
      if (ang_mom>b->max_ang_mom)
        b->max_ang_mom = ang_mom;
      b->n_prim_shells += s->n_prim;
      s->n_bfce = basis_shell_nbfce(s->type);
      c->n_bfce += s->n_bfce;
      bf_id += s->n_bfce;
      }
    b->n_bfce += c->n_bfce;
    }
  if (!b->n_ibfce)
    b->n_ibfce = b->n_bfce;
  }

/* -------------------------------------------------------------------------- */

/* define specified basis set
 
   b - pointer to basis set struct
   d - user-define basis set
   c - basis set type identifier */
void basis_def(struct basis *b, struct basis *d, short c) {
  switch (c) {
    case BASIS_DEF_STO3G:
      basis_def_sto3g(b);
      break;
    case BASIS_DEF_GEN: 
      basis_def_gen(b,d);
      break;
    }
  basis_def_parm(b);
  }

/* -------------------------------------------------------------------------- */
