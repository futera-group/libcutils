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
#include <cmn/file.h>
#include <cmn/vector.h>
#include <cmn/string.h>
#include <mol/atom.h>
#include "qmc/basis.h"

/* -------------------------------------------------------------------------- */

/* write basis set in molpro format to file
 
   b    - pointer to basis set struct
   name - name of the file */
void basis_molpro_write(struct basis *b, char *name) {
  unsigned ia,ic,is,ip,it,at[120],n_cont;
  short shell_new;
  double *cf;
  char *type;
  FILE *file = stdout;
  struct basis_center *c;
  struct basis_shell *s;
  if (!b)
    return;
  if (name && name[0])
    file = file_open(name,"w");
  /* select different atoms */
  vec_uset(at,0,120);
  for (ic=0; ic<b->n_centers; ic++)
    at[b->center[ic].type]=1;
  fprintf(file,"basis\n");
  /* one unique atom type */
  for (ia=0; ia<120; ia++)
    if (at[ia]) {
      ic = 0;
      while (b->center[ic].type!=ia)
        ic++;
      c = b->center+ic;
      /* sort by shell types */
      for (it=1; it<BASIS_SHELL_NUM; it++) {
        shell_new = 1;
        /* exponents */
        for (is=0; is<c->n_shells; is++) {
          s = c->shell+is;
          if (it!=BASIS_SHELL_SP && (s->type==it ||
             (s->type==BASIS_SHELL_SP && it==BASIS_SHELL_S) ||
             (s->type==BASIS_SHELL_SP && it==BASIS_SHELL_P))) {
            /* shell type label */
            if (s->type==BASIS_SHELL_SP && it==BASIS_SHELL_S)
              type = basis_shell_name(BASIS_SHELL_S);
            else if (s->type==BASIS_SHELL_SP && it==BASIS_SHELL_P)
              type = basis_shell_name(BASIS_SHELL_P);
            else
              type = basis_shell_name(s->type);
            str_lowcase(type);
            /* print shell type, atom name and exponents */
            if (shell_new) {
              fprintf(file,"%s, %s",type,atom_name(c->type));
              shell_new = 0;
              }
            for (ip=0; ip<s->n_prim; ip++)
              fprintf(file,", %14.10f",s->exp[ip]);
            }
          }
        if (!shell_new)
          fprintf(file,"\n");
        /* contraction coefficients */
        n_cont = 1;
        for (is=0; is<c->n_shells; is++) {
          s = c->shell+is;
          if (it!=BASIS_SHELL_SP && (s->type==it ||
             (s->type==BASIS_SHELL_SP && it==BASIS_SHELL_S) ||
             (s->type==BASIS_SHELL_SP && it==BASIS_SHELL_P))) {
            fprintf(file,"c, %d.%d",n_cont,n_cont+s->n_prim-1);
            cf = ((s->type==BASIS_SHELL_SP && it==BASIS_SHELL_P) ? 
              s->cf2 : s->cf1);
            for (ip=0; ip<s->n_prim; ip++)
              fprintf(file,", %14.10f",cf[ip]);
            fprintf(file,"\n");
            n_cont += s->n_prim;
            }
          }
        }
      }
  fprintf(file,"end\n");
  /* finish */
  if (name && name[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
