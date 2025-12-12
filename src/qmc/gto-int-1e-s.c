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
#include <cmn/math.h>
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* calculate overlap integral between two core GTO functions

   e1,e2 - GTO exponents
   n1,n2 - GTO specifiers (l,m,n)
   pa,pb - radius vectors from centers */
double gto_int_1e_s_core(double e1, unsigned *n1, double *pa,
  double e2, unsigned *n2, double *pb) {
  static double ff,sx,core;
  static unsigned i,ix;
  core = 1.0;
  for (ix=0; ix<3; ix++) {
    sx = 0.0;
    for (i=0; i<=((n1[ix]+n2[ix])/2); i++) {
      ff = gto_afce_f(2*i,n1[ix],n2[ix],pa[ix],pb[ix]);
      if (i) {
        ff *= math_fact_odd_f(2*i-1);
        ff /= pow(2.0*(e1+e2),i);
        }
      sx += ff;
      }
    core *= sx;
    }
  return(core);
  }

/* calculate overlap integral between two gaussian basis functions

   s1,s1 - pointer to basis shell GTO structs
   c1,c2 - shell center coordinates
   f1,f2 - basis functions IDs
   rr    - inter-function vector
   r2    - square of the rr vector */
double gto_int_1e_s_bf(struct basis_shell *s1, double *c1, unsigned f1,
  struct basis_shell *s2, double *c2, unsigned f2, double *rr, double r2) {
  static unsigned f1_id[3*GTO_LC_MAX],f2_id[3*GTO_LC_MAX],f1_n,f2_n;
  static double pp[3],p1[3],p2[3],nm1,nm2,cf1,cf2,pf1,pf2,pref,core;
  static double ss,ee,e12,f1_mf[GTO_LC_MAX],f2_mf[GTO_LC_MAX];
  static unsigned i,if1,if2,ip1,ip2;
  /* pure to cartesian transformation */
  gto_pure_to_cart(s1->type,f1,f1_id,f1_mf,&f1_n);
  gto_pure_to_cart(s2->type,f2,f2_id,f2_mf,&f2_n);
  /* evaluate overlap */
  ss = 0.0;
  /* loop over contracted functions */
  for (if1=0; if1<f1_n; if1++) {
    for (if2=0; if2<f2_n; if2++) {
    /* loop over primitive functions */
      for (ip1=0; ip1<s1->n_prim; ip1++) {
        for (ip2=0; ip2<s2->n_prim; ip2++) {
          /* prefactor */
          e12 = s1->exp[ip1]+s2->exp[ip2];
          for (i=0; i<3; i++) {
            pp[i] = (s1->exp[ip1]*c1[i]+s2->exp[ip2]*c2[i])/e12;
            p1[i] = pp[i]-c1[i];
            p2[i] = pp[i]-c2[i];
            }
          pf1 = pow(M_PI/e12,1.5);
          ee = s1->exp[ip1]*s2->exp[ip2]/e12;
          pf2 = exp(-ee*r2);
          /* core integral */
          core = gto_int_1e_s_core(
            s1->exp[ip1],f1_id+3*if1,p1,
            s2->exp[ip2],f2_id+3*if2,p2);
          /* normalization */
          nm1 = gto_norm(s1->exp[ip1],f1_id+3*if1);
          nm2 = gto_norm(s2->exp[ip2],f2_id+3*if2);
          /* contraction */
          cf1 = (s1->type==BASIS_SHELL_SP && f1 ? s1->cf2[ip1] : s1->cf1[ip1]);
          cf2 = (s2->type==BASIS_SHELL_SP && f2 ? s2->cf2[ip2] : s2->cf1[ip2]);
          pref = nm1*nm2*cf1*cf2*pf1*pf2*f1_mf[if1]*f2_mf[if2];
          ss += (pref*core);
          }
        }
      }
    }
  return(ss);
  }

/* calculate overlap matrix 

   b - pointer to basis shell GTO structs
   m - storage for the matrix */
void gto_int_1e_s_mat(struct basis *b, double **m) {
  unsigned i,ic1,ic2,is1,is2,if1,if2,id1,id2;
  double rr[3],r2;
  struct basis_center *c1,*c2;
  struct basis_shell *s1,*s2;
  /* loop over centers */
  for (ic1=0; ic1<b->n_centers; ic1++) {
    c1 = b->center+ic1;
    for (ic2=0; ic2<b->n_centers; ic2++) {
      c2 = b->center+ic2;
      /* loop over shells */
      for (is1=0; is1<c1->n_shells; is1++) {
        s1 = c1->shell+is1;
        for (is2=0; is2<c2->n_shells; is2++) {
          s2 = c2->shell+is2;
          /* inter-shell distance */
          r2 = 0.0;
          for (i=0; i<3; i++) {
            rr[i] = c1->coord[i]-c2->coord[i];
            r2 += rr[i]*rr[i];
            }
          /* loop over basis functions */
          for (if1=0; if1<s1->n_bfce; if1++) {
            id1 = s1->bf+if1;
            for (if2=0; if2<s2->n_bfce; if2++) {
              id2 = s2->bf+if2;
              /* matrix element */
              if (id1==id2)
                m[id1][id2] = 1.0;
              else if (id1<id2)
                m[id1][id2] = m[id2][id1]=
                 gto_int_1e_s_bf(s1,c1->coord,if1,s2,c2->coord,if2,rr,r2);
              }
            }
          }
        }
      }
    }
  }

/* calculate overlap matrix between two basis sets

   b1,b2 - pointer to basis shell GTO structs
   m     - storage for the matrix */
void gto_int_1e_s_mat_b12(struct basis *b1, struct basis *b2, double **m) {
  unsigned i,ic1,ic2,is1,is2,if1,if2,id1,id2;
  double rr[3],r2;
  struct basis_center *c1,*c2;
  struct basis_shell *s1,*s2;
  /* loop over centers */
  for (ic1=0; ic1<b1->n_centers; ic1++) {
    c1 = b1->center+ic1;
    for (ic2=0; ic2<b2->n_centers; ic2++) {
      c2 = b2->center+ic2;
      /* loop over shells */
      for (is1=0; is1<c1->n_shells; is1++) {
        s1 = c1->shell+is1;
        for (is2=0; is2<c2->n_shells; is2++) {
          s2 = c2->shell+is2;
          /* inter-shell distance */
          r2 = 0.0;
          for (i=0; i<3; i++) {
            rr[i] = c1->coord[i]-c2->coord[i];
            r2 += rr[i]*rr[i];
            }
          /* loop over basis functions */
          for (if1=0; if1<s1->n_bfce; if1++) {
            id1 = s1->bf+if1;
            for (if2=0; if2<s2->n_bfce; if2++) {
              id2 = s2->bf+if2;
              /* matrix element */
              m[id1][id2]=
                gto_int_1e_s_bf(s1,c1->coord,if1,s2,c2->coord,if2,rr,r2);
              }
            }
          }
        }
      }
    }
  }

/* -------------------------------------------------------------------------- */
