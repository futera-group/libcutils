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
#include <cmn/string.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include "qmc/basis.h"

/* -------------------------------------------------------------------------- */

/* write basis set in turbomole format to file 
 
   b    - pointer to basis set struct
   name - name of the file */
void basis_turbo_write(struct basis *b, char *name) {
  unsigned ia,ic,is,ip,at[120];
  char *sym,type[20] = "\0";
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
    at[b->center[ic].type] = 1;
  fprintf(file,"$basis\n*\n");
  for (ia=0; ia<120; ia++)
    if (at[ia]) {
      sym = atom_name(ia);
      str_lowcase(sym);
      /* print atom name */
      fprintf(file,"%s BasisSetName\n*\n",sym);
      ic=0;
      while (b->center[ic].type!=ia)
        ic++;
      c = b->center+ic;
      for (is=0; is<c->n_shells; is++) {
        s = c->shell+is;
        if (s->type==BASIS_SHELL_SP) 
          sprintf(type,"s");
        else {
          sprintf(type,"%s",basis_shell_name(s->type));
          str_lowcase(type);
          }
        fprintf(file,"%3d  %s\n",s->n_prim,type);
        for (ip=0; ip<s->n_prim; ip++) 
          fprintf(file,"%18.10f %18.10f\n",s->exp[ip],s->cf1[ip]);
        if (s->type==BASIS_SHELL_SP) {
          fprintf(file,"%3d  p\n",s->n_prim);
          for (ip=0; ip<s->n_prim; ip++) 
            fprintf(file,"%18.10f %18.10f\n",s->exp[ip],s->cf2[ip]);
          }
        }
      fprintf(file,"*\n");
      }
  fprintf(file,"$end\n");
  /* finish */
  if (name && name[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
