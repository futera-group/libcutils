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
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* calculate one element of HF 2e integral matrix
 
   b         - pointer to basis shell GTO structs
   p         - density matrix
   s1,s2     - basis set shells
   c1,c1     - basis set shell coordinates
   f1,f2     - basis function IDs
   g_ij,g_ji - G matrix elements
   up        - criterion for neglecting pair contributions
   uu        - criterion for neglecting pair-pair contributions
   zr        - criterion for neglecting the matrix element */
void gto_int_2e_g_elm(struct basis *b, double **p,
  struct basis_shell *s1, double *c1, unsigned f1,
  struct basis_shell *s2, double *c2, unsigned f2,
  double *g_ij, double *g_ji, double up, double uu, double zr) {
  static unsigned ic3,ic4,is3,is4,if3,if4,id1,id2,id3,id4;
  static double t1,t2,t3,tt;
  static struct basis_center *c3,*c4;
  static struct basis_shell *s3,*s4;
  id1 = s1->bf+f1;
  id2 = s2->bf+f2;
  (*g_ij) = (*g_ji)=0.0;
  /* third center */
  for (ic3=0; ic3<b->n_centers; ic3++) {
    c3 = b->center+ic3;
    for (is3=0; is3<c3->n_shells; is3++) {
      s3 = c3->shell+is3;
      for (if3=0; if3<s3->n_bfce; if3++) {
        id3 = s3->bf+if3;
        /* fourth center */
        for (ic4=0; ic4<b->n_centers; ic4++) {
          c4 = b->center+ic4;
          for (is4=0; is4<c4->n_shells; is4++) {
            s4 = c4->shell+is4;
            for (if4=0; if4<s4->n_bfce; if4++) {
              id4 = s4->bf+if4;
              /* matrix element */
              if (fabs(p[id3][id4])>zr) {
                if (id3>=id4) {
                  t1 = gto_int_2e_e_bf_r(
                    s1,c1,f1,s2,c2,f2,
                    s3,c3->coord,if3,s4,c4->coord,if4,up,uu,zr);
                  t2 = gto_int_2e_e_bf_r(
                    s1,c1,f1,s3,c3->coord,if3,
                    s4,c4->coord,if4,s2,c2,f2,up,uu,zr);
                  if (id3==id4) {
                    tt = (t1-0.5*t2);
                    (*g_ij) += (p[id3][id4]*tt);
                    if (id1!=id2)
                      (*g_ji) += (p[id4][id3]*tt);
                    }
                  else {
                    tt = (2.0*t1-0.5*t2);
                    if (id1==id2)
                      (*g_ij) += (p[id3][id4]*tt-0.5*p[id4][id3]*t2);
                    else {
                      t3 = gto_int_2e_e_bf_r(
                        s1,c1,f1,s4,c4->coord,if4,
                        s3,c3->coord,if3,s2,c2,f2,up,uu,zr);
                      (*g_ij) += (p[id3][id4]*tt-0.5*p[id4][id3]*t3);
                      (*g_ji) += (p[id4][id3]*tt-0.5*p[id3][id4]*t3);
                      }
                    }
                  }
                else
                  break;
                }
              if (id3<id4)
                break;
              }
            if (id3<id4)
              break;
            }
          if (id3<id4)
            break;
          }
        }
      }
    }
  }

/* calculate HF 2e integral matrix
 
   b  - pointer to basis shell GTO structs
   p  - density matrix
   g  - HF 2e integral matrix
   up - criterion for neglecting pair contributions
   uu - criterion for neglecting pair-pair contributions
   zr - criterion for neglecting the matrix element */
void gto_int_2e_g_mat(struct basis *b, double **p, double **g,
  double up, double uu, double zr) {
  unsigned ic1,ic2,is1,is2,if1,if2,id1,id2=0;
  struct basis_center *c1,*c2;
  struct basis_shell *s1,*s2;
  /* first center */
  for (ic1=0; ic1<b->n_centers; ic1++) {
    c1 = b->center+ic1;
    for (is1=0; is1<c1->n_shells; is1++) {
      s1 = c1->shell+is1;
      for (if1=0; if1<s1->n_bfce; if1++) {
        id1 = s1->bf+if1;
        /* second center */
        for (ic2=0; ic2<b->n_centers; ic2++) {
          c2 = b->center+ic2;
          for (is2=0; is2<c2->n_shells; is2++) {
            s2 = c2->shell+is2;
            for (if2=0; if2<s2->n_bfce; if2++) {
              id2 = s2->bf+if2;
              /* matrix element */
              if (id1>=id2)
                gto_int_2e_g_elm(b,p,s1,c1->coord,if1,s2,c2->coord,if2,
                  &(g[id1][id2]),&(g[id2][id1]),up,uu,zr);
              else
                break;
              }
            if (id1<id2)
              break;
            }
          if (id1<id2)
            break;
          }
        }
      }
    }
  }

/* -------------------------------------------------------------------------- */
