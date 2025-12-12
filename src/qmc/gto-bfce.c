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

#include <math.h>
#include "qmc/basis.h"
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* calculate values of all basis functions at specific point

   c  - pointer to basis set center data struct
   s  - pointer to basis set shell data struct
   id - ID of the basis function
   x  - coordinates of the point */
double gto_bfce_val(struct basis_center *c, struct basis_shell *s, unsigned id,
  double *x) {
  unsigned i,*a,f_id[3*GTO_LC_MAX],f_n,ig,ip;
  double val,bf,rr,r2,f_mf[GTO_LC_MAX];
  /* pure to cartesian transformation */
  gto_pure_to_cart(s->type,id,f_id,f_mf,&f_n);
  /* function evaluation */
  bf = 0.0;
  for (ig=0; ig<f_n; ig++) {
    a = f_id+3*ig;
    for (ip=0; ip<s->n_prim; ip++) {
      val = 1.0;
      r2 = 0.0;
      for (i=0; i<3; i++) {
        rr = x[i]-c->coord[i];
        r2 += rr*rr;
        if (a[i])
          val *= pow(rr,a[i]);
        }
      val *= exp(-(s->exp[ip]*r2));
      val *= gto_norm(s->exp[ip],a);
      val *= f_mf[ig];
      val *= (s->type==BASIS_SHELL_SP && id ? s->cf2[ip] : s->cf1[ip]);
      bf += val;
      }
    }
  return(bf);
  }

/* calculate values of all basis functions at specific point

   b - pointer to basis set data struct
   x - coordinates of the point
   f - storage array for the value of basis function */
void gto_bfce_val_all(struct basis *b, double *x, double *f) {
  unsigned ic,is,ib,bf_id;
  struct basis_center *c;
  struct basis_shell *s;
  bf_id = 0;
  for (ic=0; ic<b->n_centers; ic++) {
    c = b->center+ic;
    for (is=0; is<c->n_shells; is++) {
      s = c->shell+is;
      for (ib=0; ib<s->n_bfce; ib++) {
        f[bf_id] = gto_bfce_val(c,s,ib,x);
        bf_id++;
        }
      }
    }
  }

/* -------------------------------------------------------------------------- */
