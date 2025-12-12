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
#include <cmn/message.h>
#include "qmc/basis.h"
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* calculate argument for Boys function
 
   p1,p2 - GTO pair data
   pq    - P-Q vector between bra and ket centers
   th    - combination of exponential factors (theta)
   tt    - argument for Boys function */
void gto_int_2e_e_parm(struct gto_pair *p1, struct gto_pair *p2,
  double *pq, double *th, double *tt) {
  static unsigned i;
  static double pq2;
  pq2 = 0.0;
  for (i=0; i<3; i++) {
    pq[i] = (p2->crd[i]-p1->crd[i]);
    pq2 += (pq[i]*pq[i]);
    }
  (*th) = (0.5/(p1->sigma+p2->sigma));
  (*tt) = ((*th)*pq2);
  }

/* evaluate one two electron integral from pre-calculated data
 
   p1,p2 - GTO pair data
   pq    - vector P-Q between centers of bra and ket
   uu    - product of bra and ket integral prefactor
   th    - combination of exponential factors (theta)
   tt    - argument for Boys function
   zr    - criterion for neglecting the contribution */
double gto_int_2e_e_core(struct gto_pair *p1, struct gto_pair *p2,
  double *pq, double uu, double th, double tt, double zr) {
  static double rm[2*(GTO_MAX_AMOM+1)+1][2*(GTO_MAX_AMOM+1)+1]
                  [2*(GTO_MAX_AMOM+1)+1][4*(GTO_MAX_AMOM+1)+1];
  static unsigned i,j,k,q,r[3],l_tot;
  static double s1,s2,e2,vv;
  l_tot = p1->am1+p1->am2+p2->am1+p2->am2;
  /* evaluate set of [r]^m integrals */
  if (!gto_int_2e_rm(rm,l_tot,uu,th,tt,pq))
    msg_error_f("[%d]^m integrals needed for %s are not implemented",
      1,l_tot-1,gto_int_2e_name(p1->am1,p1->am2,p2->am1,p2->am2));
  /* expand bra and ket to basic terms */
  gto_int_2e_pq(p1,zr);
  gto_int_2e_pq(p2,zr);
  /* final integration */
  e2 = 0.0;
  for (i=0; i<p1->n_pq; i++)
    for (j=0; j<p2->n_pq; j++) {
      q = 0;
      for (k=0; k<3; k++) {
        r[k] = p1->pq[i].id[k]+p2->pq[j].id[k];
        q += p2->pq[j].id[k];
        }
      if ((r[0]+r[1]+r[2])>l_tot)
        msg_error_f("r[%d]^0 required for %s but maximal r[%d]^m types"
          " are pre-calculated",1,r[0]+r[1]+r[2],
          gto_int_2e_name(p1->am1,p1->am2,p2->am1,p2->am2),l_tot);
      /* pq-braket --> (-1)^q*[r]^0 */
      vv = rm[r[0]][r[1]][r[2]][0];
      if (q%2)
        vv = (-vv);
      s1 = (pow(2.0*p1->e12,p1->am1+p1->am2)*p1->pq[i].coeff);
      s2 = (pow(2.0*p2->e12,p2->am1+p2->am2)*p2->pq[j].coeff);
      e2 += (s1*s2*vv);
      }
  return(e2);
  }

/* evaluate one (ab|cd) 2e integral from 4 contracted GTO shells
 
   s1,s2,s3,s4 - basis set shells
   c1,c2,c3,c4 - basis set shell coordinates
   f1,f2,f3,f4 - basis set function IDs
   r2_12,r2_34 - square distances between basis shell pairs
   r12,r34     - intershell vectors
   r2_12,r2_34 - square intershell distances
   up          - criterion for neglecting pair contributions
   uu          - criterion for neglecting pair-pair contributions
   zr          - criterion for neglecting the contribution */
double gto_int_2e_bf(
  struct basis_shell *s1, double *c1, unsigned f1,
  struct basis_shell *s2, double *c2, unsigned f2,
  struct basis_shell *s3, double *c3, unsigned f3,
  struct basis_shell *s4, double *c4, unsigned f4,
  double *r12, double r2_12, double *r34, double r2_34,
  double up, double uu, double zr) {
  static unsigned f1_id[3*GTO_LC_MAX],f2_id[3*GTO_LC_MAX],f1_n,f2_n;
  static unsigned f3_id[3*GTO_LC_MAX],f4_id[3*GTO_LC_MAX],f3_n,f4_n;
  static double f1_mf[GTO_LC_MAX],f2_mf[GTO_LC_MAX],integral,pq[3],th,tt;
  static double f3_mf[GTO_LC_MAX],f4_mf[GTO_LC_MAX],pref;
  static unsigned ip12,ip34,if1,if2,if3,if4;
  static long unsigned n12,n34,n12m,n34m;
  static struct gto_pair *p12,*p34;
  /* pure to cartesian transformation */
  gto_pure_to_cart(s1->type,f1,f1_id,f1_mf,&f1_n);
  gto_pure_to_cart(s2->type,f2,f2_id,f2_mf,&f2_n);
  gto_pure_to_cart(s3->type,f3,f3_id,f3_mf,&f3_n);
  gto_pure_to_cart(s4->type,f4,f4_id,f4_mf,&f4_n);
  /* arrays for pair data */
  n12m = s1->n_prim*s2->n_prim*f1_n*f2_n;
  n34m = s3->n_prim*s4->n_prim*f3_n*f4_n;
  p12 = gto_pair_new(n12m);
  p34 = gto_pair_new(n34m);
 /* pair data generation */
  n12 = n34 = 0;
  for (if1=0; if1<f1_n; if1++)
    for (if2=0; if2<f2_n; if2++)
      gto_pair_eval_one(p12,
        s1,c1,f1,f1_id+3*if1,
        s2,c2,f2,f2_id+3*if2,
        r12,r2_12,f1_mf[if1]*f2_mf[if2],
        n12,&n12,n12m,up);
  for (if3=0; if3<f3_n; if3++)
    for (if4=0; if4<f4_n; if4++)
      gto_pair_eval_one(p34,
        s3,c3,f3,f3_id+3*if3,
        s4,c4,f4,f4_id+3*if4,
        r34,r2_34,f3_mf[if3]*f4_mf[if4],
        n34,&n34,n34m,up);
  /* select integrals and make sum of them */
  integral = 0.0;
  for (ip12=0; ip12<n12; ip12++)
    for (ip34=0; ip34<n34; ip34++) {
      pref = p12[ip12].pref*p34[ip34].pref;
      if (fabs(pref)>uu) {
        /* set parameters of Boys function */
        gto_int_2e_e_parm(&(p12[ip12]),&(p34[ip34]),pq,&th,&tt);
        /* (ab|cd) integral */
        integral += 
          gto_int_2e_e_core(&(p12[ip12]),&(p34[ip34]),pq,pref,th,tt,zr);
        }
      }
  /* clean memory */
  gto_pair_free(p12);
  gto_pair_free(p34);
  return(integral);
  }

/* evaluate one (ab|cd) 2e integral from 4 contracted GTO shells
 
   s1,s2,s3,s4 - basis set shells
   c1,c2,c3,c4 - basis set shell coordinates
   f1,f2,f3,f4 - basis set function IDs
   up          - criterion for neglecting pair contributions
   uu          - criterion for neglecting pair-pair contributions
   zr          - criterion for neglecting the contribution */
double gto_int_2e_e_bf_r(
  struct basis_shell *s1, double *c1, unsigned f1,
  struct basis_shell *s2, double *c2, unsigned f2,
  struct basis_shell *s3, double *c3, unsigned f3,
  struct basis_shell *s4, double *c4, unsigned f4,
  double up, double uu, double zr) {
  static double r12[3],r34[3],r2_12,r2_34;
  static unsigned i;
  r2_12 = r2_34 = 0.0;
  for (i=0; i<3; i++) {
    r12[i] = c1[i]-c2[i];
    r2_12 += r12[i]*r12[i];
    r34[i] = c3[i]-c4[i];
    r2_34 += r34[i]*r34[i];
    }
  return(gto_int_2e_bf(
    s1,c1,f1,s2,c2,f2,
    s3,c3,f3,s4,c4,f4,
    r12,r2_12,r34,r2_34,up,uu,zr));
  }

/* -------------------------------------------------------------------------- */
