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
#include <cmn/message.h>
#include "qmc/basis.h"
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* create copy of GTO pair with swapped data
 
   a,b - pointers to the GTO pair data structs */
void gto_pair_swap(struct gto_pair *a, struct gto_pair *b) {
  unsigned i;
  /* exponents */
  a->e1 = b->e2;
  a->e2 = b->e1;
  a->e12 = b->e12;
  /* angular momenta */
  a->am1 = b->am2;
  a->am2 = b->am1;
  /* vectors */
  for (i=0; i<3; i++) {
    a->aa1[i] = b->aa2[i];
    a->aa2[i] = b->aa1[i];
    a->crd[i] = b->crd[i];
    a->r12[i] = (-b->r12[i]);
    }
  a->sigma = b->sigma;
  a->pref = b->pref;
  a->n_pq = 0;
  }

/* -------------------------------------------------------------------------- */

/* evaluate core data for one GTO basis function pair

   p     - GTO pair data struct array
   a1,a2 - GTO angular momenta (modulus)
   c1,c2 - GTO coordinates
   e1,e2 - GTO exponents 
   d1,d2 - GTO contraction coefficients
   r2    - square of intershell distance */
void gto_pair_eval_core(struct gto_pair *p,
  unsigned a1, double *c1, double e1, double d1, 
  unsigned a2, double *c2, double e2, double d2,
  double *r12, double r2) {
  static unsigned i;
  /* exponents */
  p->e1 = e1;
  p->e2 = e2;
  p->e12 = e1+e2;
  /* angular momenta */
  p->am1 = a1;
  p->am2 = a2;
  /* sigma */
  p->sigma = 1.0/(2.0*p->e12);
  /* pair center */
  for (i=0; i<3; i++) {
    p->r12[i] = r12[i];
    p->crd[i] = (e1*c1[i]+e2*c2[i])/p->e12;
    }
  /* prefactor */
  p->pref = MATH_2r8PIs3;
  p->pref *= pow(p->sigma,p->am1+p->am2+1.5);
  p->pref *= (d1*d2);
  p->pref *= exp(-e1*e2*r2/p->e12);
  }

/* evaluated one pair GTO data
 
   p     - GTO pair data struct array
   s1,s2 - basis set shells
   c1,c2 - basis set shell coordinates
   f1,f2 - basis set function IDs
   a1,a2 - angular momenta of the functions
   r2    - square of intershell distance
   n0,n1 - initial and final number of pairs
   np    - maximum number of pairs
   zero  - criterion for neglecting the pair */
void gto_pair_eval_one(struct gto_pair *p,
  struct basis_shell *s1, double *c1, unsigned f1, unsigned *a1,
  struct basis_shell *s2, double *c2, unsigned f2, unsigned *a2,
  double *r12, double r2, double sc, long unsigned n0,
  long unsigned *n1, long unsigned np, double zero) {
  static unsigned i,ip1,ip2,a1loc[3],a2loc[3];
  static long unsigned id;
  id = n0;
  /* angular momentum vectors */
  if (!a1) {
    gto_ang_vec(s1->type,f1,a1loc);
    a1 = a1loc;
    }
  if (!a2) {
    gto_ang_vec(s2->type,f2,a2loc);
    a2 = a2loc;
    }
  /* loop over primitive functions */
  for (ip1=0; ip1<s1->n_prim; ip1++) {
    for (ip2=0; ip2<s2->n_prim; ip2++) {
      if (id>=np)
        msg_error_f("maximum number of GTO pairs exceeded (%d)",1,np);
      /* core pair data */
      gto_pair_eval_core(&(p[id]),
        a1[0]+a1[1]+a1[2],c1,s1->exp[ip1],
        (s1->type==BASIS_SHELL_SP && f1 ? s1->cf2[ip1] : s1->cf1[ip1]),
        a2[0]+a2[1]+a2[2],c2,s2->exp[ip2],
        (s2->type==BASIS_SHELL_SP && f2 ? s2->cf2[ip2] : s2->cf1[ip2]),r12,r2);
      /* normalization constants */
      p[id].pref *= gto_norm(s1->exp[ip1],a1);
      p[id].pref *= gto_norm(s2->exp[ip2],a2);
      /* pre-selection */
      p[id].pref *= sc;
      if (fabs(p[id].pref)>=zero) {
        for (i=0; i<3; i++) {
          p[id].aa1[i] = a1[i];
          p[id].aa2[i] = a2[i];
          }
        id++;
        }
      }
    }
  (*n1) = id;
  }

/* calculate and GTO pair integral data and select only important terms
 
   ab    - pointer to GTO pair data array 
   b1,b2 - basis sets
   m1,m2 - MO coefficients
   np    - number of selected pairs (output)
   nt    - total number of allocated pairs
   up    - criterion for neglecting pair contributions */
void gto_pair_eval_all(struct gto_pair *ab, struct basis *b1, double *m1,
  struct basis *b2, double *m2, long unsigned *np, long unsigned nt,
  double up) {
  unsigned i,ic1,ic2,is1,is2,if1,if2,id1,id2,ip1,ip2,*a1,*a2;
  unsigned f1_id[3*GTO_LC_MAX],f2_id[3*GTO_LC_MAX],f1_n,f2_n;
  double f1_mf[GTO_LC_MAX],f2_mf[GTO_LC_MAX],scl,r12[3],r2;
  struct basis_center *c1,*c2;
  struct basis_shell *s1,*s2;
  long unsigned n0,n1;
  n0 = 0;
  /* loop over centers */
  for (ic1=0; ic1<b1->n_centers; ic1++) {
    c1 = b1->center+ic1;
    for (ic2=0; ic2<b2->n_centers; ic2++) {
      c2 = b2->center+ic2;
      /* intershell distance */
      r2 = 0.0;
      for (i=0; i<3; i++) {
        r12[i] = c1->coord[i]-c2->coord[i];
        r2 += (r12[i]*r12[i]);
        }
      /* loop over shells */
      for (is1=0; is1<c1->n_shells; is1++) {
        s1 = c1->shell+is1;
        for (is2=0; is2<c2->n_shells; is2++) {
          s2 = c2->shell+is2;
          /* loop over basis functions */
          for (if1=0; if1<s1->n_bfce; if1++) {
            id1 = s1->bf+if1;
            gto_pure_to_cart(s1->type,if1,f1_id,f1_mf,&f1_n);
            for (if2=0; if2<s2->n_bfce; if2++) {
              id2 = s2->bf+if2;
              gto_pure_to_cart(s2->type,if2,f2_id,f2_mf,&f2_n);
              /* Talmi transformation */
              for (ip1=0; ip1<f1_n; ip1++) {
                a1 = f1_id+3*ip1;
                for (ip2=0; ip2<f2_n; ip2++) {
                  a2 = f2_id+3*ip2;
                  /* pair data calculation */
                  scl = (m1[id1]*m2[id2]);
                  scl *= (f1_mf[ip1]*f2_mf[ip2]);
                  gto_pair_eval_one(ab,
                    s1,c1->coord,if1,a1,
                    s2,c2->coord,if2,a2,
                    r12,r2,scl,n0,&n1,nt,up);
                  n0 = n1;
                  }
                }
              }
            }
          }
        }
      }
    }
  (*np) = n1;
  }

/* -------------------------------------------------------------------------- */
