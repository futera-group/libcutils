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
#include <mol/atom.h>
#include "qmc/basis.h"

/* -------------------------------------------------------------------------- */

/* print basis set info to file
 
   b    - pointer to basis set struct
   t    - type of basis set print-out
   file - pointer to open output stream */
void basis_fprint(struct basis *b, short t, FILE *file) {
  unsigned ic,is,ip,id_s;
  struct basis_center *c;
  struct basis_shell *s;
  if (!b)
    return;
  if (t==2) {
    /* centers */
    id_s = 1;
    for (ic=0; ic<b->n_centers; ic++) {
      c = b->center+ic;
      fprintf(file,"Center %3d (%-3s) %12.6f %12.6f %12.6f\n",
        ic+1,atom_name(c->type),c->coord[0],c->coord[1],c->coord[2]);
      for (is=0; is<c->n_shells; is++) {
        s = c->shell+is;
        fprintf(file,"  Shell %4d (%-2s)  NPrim %2d  BFce %5d - %5d\n",
          id_s,basis_shell_name(s->type),s->n_prim,s->bf+1,s->bf+s->n_bfce);
        for (ip=0; ip<s->n_prim; ip++) {
          fprintf(file,"%20.10E %18.10E",s->exp[ip],s->cf1[ip]);
          if (s->type==BASIS_SHELL_SP)
            fprintf(file," %18.10E",s->cf2[ip]);
          fprintf(file,"\n");
          }
        id_s++;
        }
      }
    }
  else {
    /* shells */
    id_s = 1;
    for (ic=0; ic<b->n_centers; ic++) {
      c = b->center+ic;
       for (is=0; is<c->n_shells; is++) {
         s = c->shell+is;
         fprintf(file,"Center %3d (%-3s) Shell %4d (%-2s) %2d %12.6f %12.6f" 
           " %12.6f\n",ic+1,atom_name(c->type),id_s,
           basis_shell_name(s->type),s->n_prim,
           c->coord[0],c->coord[1],c->coord[2]);
         for (ip=0; ip<s->n_prim; ip++) {
           fprintf(file,"%18.10E %18.10E",s->exp[ip],s->cf1[ip]);
           if (s->type==BASIS_SHELL_SP)
             fprintf(file," %18.10E",s->cf2[ip]);
           fprintf(file,"\n");
          }
        id_s++;
        }
      }
    }
  }

/* print basis set info 
 
   b - pointer to basis set struct
   t - type of basis set print-out
   f - name of the output file */
void basis_print(struct basis *b, short t, char *f) {
  FILE *file = stdout;
  if (!b)
    return;
  if (f && f[0])
    file = file_open(f,"w");
  basis_fprint(b,t,file);
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
