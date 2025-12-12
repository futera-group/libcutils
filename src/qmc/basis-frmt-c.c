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

#define BASIS_PRIM_EXP  1  /* primitive shell exponents */
#define BASIS_PRIM_CF1  2  /* primitive shell coefficients */
#define BASIS_PRIM_CF2  3  /* coefficients of P functions in SP shells */

/* -------------------------------------------------------------------------- */

/* write one basis set array in C format to file 
 
   b      - pointer to basis set struct
   at     - array with atom indicators 
   n_prim - total number of primitive shells
   file   - pointer to open output stream
   type   - type of output (exp, cf1 or cf2) */
void basis_c_write_one(struct basis *b, unsigned *at, unsigned n_prim,
  FILE *file, short type) {
  unsigned ia,ic,is,ip,n_tot = 0;
  struct basis_center *c;
  struct basis_shell *s;
  for (ia=0; ia<120; ia++)
    if (at[ia]) {
      ic = 0;
      while (b->center[ic].type!=ia)
        ic++;
      c = b->center+ic;
      for (is=0; is<c->n_shells; is++) {
        s = c->shell+is;
        n_tot += b->center[ic].shell[is].n_prim;
        for (ip=0; ip<s->n_prim; ip++) {
          fprintf(file,(ip ? "%17.10E" : "%18.10E"),
            (type==BASIS_PRIM_EXP ? s->exp[ip] :
            (type==BASIS_PRIM_CF1 ? s->cf1[ip] :
            (s->type==BASIS_SHELL_SP ? s->cf2[ip] : 0.0))));
          if (ip!=(s->n_prim-1))
            fprintf(file,",");
          }
        fprintf(file,(n_tot<n_prim ? ",  " : "}; "));
        fprintf(file,"/* %-2s - %-2s */",atom_name(ia),
          basis_shell_name(s->type));
        fprintf(file,"\n");
        }
      }
  }

/* write basis set in C format to file 
 
   b    - pointer to basis set struct
   name - name of the file */
void basis_c_write(struct basis *b, char *name) {
  unsigned ia,ic,is,at[120],n_atom = 0,n_prim = 0,id_a,id_i;
  FILE *file = stdout;
  struct basis_center *c;
  if (!b)
    return;
  if (name && name[0])
    file = file_open(name,"w");
  /* select different atoms */
  vec_uset(at,0,120);
  for (ic=0; ic<b->n_centers; ic++) {
    if (!at[b->center[ic].type]) {
      for (is=0; is<b->center[ic].n_shells; is++)
        n_prim+=b->center[ic].shell[is].n_prim;
      n_atom++;
      }
    at[b->center[ic].type] = 1;
    }
  /* print out array of indices */
  id_a = 0;
  id_i = 1;
  fprintf(file,"double basis_sto3g_id[%d]={\n",n_atom);
  for (ia=0; ia<120; ia++)
    if (at[ia]) {
      ic = 0;
      while (b->center[ic].type!=ia)
        ic++;
      c = b->center+ic;
      fprintf(file,"%4d",id_a);
      for (is=0; is<c->n_shells; is++)
        id_a += c->shell[is].n_prim;
      fprintf(file,(id_a<n_prim ? "," : "};"));
      if (!(id_i%10))
        fprintf(file,"\n");
      id_i++;
      }
  if (id_i && (id_i-1)%10)
    fprintf(file,"\n");
  fprintf(file,"\n");
  /* print out array of exponents */
  fprintf(file,"double basis_sto3g_exp[%d]={\n",n_prim);
  basis_c_write_one(b,at,n_prim,file,BASIS_PRIM_EXP);
  /* print out array of coefficients */
  fprintf(file,"\n");
  fprintf(file,"double basis_sto3g_cf1[%d]={\n",n_prim);
  basis_c_write_one(b,at,n_prim,file,BASIS_PRIM_CF1);
  /* print out array of SP coefficients */
  fprintf(file,"\n");
  fprintf(file,"double basis_sto3g_cf2[%d]={\n",n_prim);
  basis_c_write_one(b,at,n_prim,file,BASIS_PRIM_CF2);
  /* finish */
  if (name && name[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
