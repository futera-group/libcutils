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
#include <cmn/matrix.h>
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* calculate potential energy integral between two core GTO functions

   n1,n2 - GTO specifiers (l,m,n)
   pa,pb - radius vectors from centers
   e12   - sum of GTO exponents
   cp    - atom center - basis sheel center vector
   cp2   - square of the cp vector */
double gto_int_1e_v_core(unsigned *n1, double *pa,
  unsigned *n2, double *pb, double e12, double *cp, double cp2) {
  static double fs1,fs2,fs3,fv1,fv2,fv3,fb;
  static unsigned i,j,k,r,s,t,u,v,w,n;
  fs1 = 0.0;
  for (i=0; i<=(n1[0]+n2[0]); i++) {
    for (r=0; r<=(i/2); r++) {
      for (u=0; u<=((i-2*r)/2); u++) {
        fv1 = gto_afce_a(i,r,u,n1[0],n2[0],pa[0],pb[0],cp[0],e12);
        fs2 = 0.0;
        for (j=0; j<=(n1[1]+n2[1]); j++) {
          for (s=0; s<=(j/2); s++) {
            for (v=0; v<=((j-2*s)/2); v++) {
              fv2 = gto_afce_a(j,s,v,n1[1],n2[1],pa[1],pb[1],cp[1],e12);
              fs3 = 0.0;
              for (k=0; k<=(n1[2]+n2[2]); k++) {
                for (t=0; t<=(k/2); t++) {
                  for (w=0; w<=((k-2*t)/2); w++) {
                    fv3 = gto_afce_a(k,t,w,n1[2],n2[2],pa[2],pb[2],cp[2],e12);
                    n = i+j+k-2*(r+s+t)-u-v-w;
                    fb = math_boys(n,e12*cp2);
                    fs3 += (fv3*fb);
                    }
                  }
                }
              fs2 += (fv2*fs3);
              }
            }
          }
        fs1 += (fv1*fs2);
        }
      }
    }
  return(fs1);
  }

/* calculate potential energy integral between two gaussian basis functions

   s1,s1 - pointer to basis shell GTO structs
   c1,c2 - shell center coordinates
   f1,f2 - basis functions IDs
   rr    - inter-function vector
   r2    - square of the rr vector
   ac    - atomic coorinates
   ah    - atomic charge */
double gto_int_1e_v_bf(struct basis_shell *s1, double *c1, unsigned f1,
  struct basis_shell *s2, double *c2, unsigned f2, double *rr, double r2,
  double *ac, double ah) {
  static double pp[3],p1[3],p2[3],cp[3],nm1,nm2,cf1,cf2,pf1,pf2,cp2,pref,core;
  static unsigned f1_id[3*GTO_LC_MAX],f2_id[3*GTO_LC_MAX],f1_n,f2_n;
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
          cp2 = 0.0;
          e12 = s1->exp[ip1]+s2->exp[ip2];
          for (i=0; i<3; i++) {
            pp[i] = (s1->exp[ip1]*c1[i]+s2->exp[ip2]*c2[i])/e12;
            p1[i] = pp[i]-c1[i];
            p2[i] = pp[i]-c2[i];
            cp[i] = pp[i]-ac[i];
            cp2 += cp[i]*cp[i];
            }
          pf1 = (2.0*M_PI)/e12;
          ee = s1->exp[ip1]*s2->exp[ip2]/e12;
          pf2 = exp(-ee*r2);
          /* core integral */
          core = gto_int_1e_v_core(f1_id+3*if1,p1,f2_id+3*if2,p2,e12,cp,cp2);
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
  return(-ah*ss);
  }

/* calculate potential energy integral between all gaussian basis functions

   b - pointer to basis set data struct 
   x - point in space where the integrals are evaluated
   v - storage array for the integrals */
void gto_int_1e_v_bf_all(struct basis *b, double *x, double **v) {
  unsigned i,j,ic1,ic2,is1,is2,if1,if2,id1,id2;
  double rr[3],r2;
  struct basis_center *c1,*c2;
  struct basis_shell *s1,*s2;
  /* initialization */
  for (i=0; i<b->n_ibfce; i++)
    for (j=0; j<b->n_ibfce; j++)
      v[i][j] = 0.0;
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
              if (id1<=id2) {
                v[id1][id2] += 
                 gto_int_1e_v_bf(s1,c1->coord,if1,s2,c2->coord,if2,
                   rr,r2,x,1.0);
                if (id1!=id2)
                  v[id2][id1] = v[id1][id2];
                }
              }
            }
          }
        }
      }
    }
  }

/* calculate potential energy matrix 

   b - pointer to basis shell GTO structs
   c - array with atomic coordinates
   h - array with atomic charges
   n - number of atoms
   m - shell center coordinates */
void gto_int_1e_v_mat(struct basis *b, double **c, double *h, unsigned n,
  double **m) {
  unsigned i,ia,ic1,ic2,is1,is2,if1,if2,id1,id2;
  double rr[3],r2;
  struct basis_center *c1,*c2;
  struct basis_shell *s1,*s2;
  /* initialization */
  mat_fset(m,0.0,b->n_ibfce,b->n_ibfce);
  /* loop over atoms */
  for (ia=0; ia<n; ia++) {
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
                if (id1<=id2) {
                  m[id1][id2] += 
                   gto_int_1e_v_bf(s1,c1->coord,if1,s2,c2->coord,if2,
                     rr,r2,c[ia],h[ia]);
                  if (id1!=id2)
                    m[id2][id1] = m[id1][id2];
                  }
                }
              }
            }
          }
        }
      }
    }
  }

/* -------------------------------------------------------------------------- */
