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
#include <mol/atom.h>
#include "qmc/basis.h"

/* -------------------------------------------------------------------------- */

/* write basis set in gamess format to file
 
   b    - pointer to basis set struct
   name - name of the file */
void basis_gamess_write(struct basis *b, char *name) {
  unsigned ia,ic,is,ip,at[120];
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
  for (ia=0; ia<120; ia++)
    if (at[ia]) {
      /* print atom name */
      fprintf(file,"%s BasisSetName\n",atom_name(ia));
      ic = 0;
      while (b->center[ic].type!=ia)
        ic++;
      c = b->center+ic;
      for (is=0; is<c->n_shells; is++) {
        s = c->shell+is;
        /* print shell type and number of primitives */
        fprintf(file,"%-2s %3d\n",
          (s->type==BASIS_SHELL_SP ? "L" : basis_shell_name(s->type)),
          s->n_prim);
        for (ip=0; ip<s->n_prim; ip++) {
          /* print coefficients and exponents of primitives */
          fprintf(file,"%3d %18.10f %18.10f",ip+1,s->exp[ip],s->cf1[ip]);
          if (s->type==BASIS_SHELL_SP)
            fprintf(file," %18.10f",s->cf2[ip]);
          fprintf(file,"\n");
          }
        }
      fprintf(file,"\n");
      }
  /* finish */
  if (name && name[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
