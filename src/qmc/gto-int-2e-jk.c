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

/* coulombic and exchange two electron integrals for one MO pair

   b1,b2 - basis sets
   m1,m2 - MO coefficients
   gj    - coulombic integral
   gk    - exchange integral
   up    - criterion for neglecting pair contributions
   pp    - criterion for neglecting pair-pair contributions
   zr    - criterion for neglecting the contribution */
void gto_calc_2e_jk(struct basis *b1, double *m1, struct basis *b2, double *m2, 
  double *gj, double *gk, double up, double pp, double zr) {
  double uu,th,tt,pq[3],ee_int,gj_tot=0.0,gk_tot=0.0;
  long unsigned n_pair;
  unsigned ip1,ip2;
  struct gto_pair ba,*ab;
  /* selection of pairs */
  ab = gto_pair_new(8*b1->n_bfce*b2->n_bfce);
  gto_pair_eval_all(ab,b1,m1,b2,m2,&n_pair,8*b1->n_bfce*b2->n_bfce,up);
  /* select integrals */
  for (ip1=0; ip1<n_pair; ip1++)
    for (ip2=0; ip2<n_pair; ip2++) {
      uu = ab[ip1].pref*ab[ip2].pref;
      if (fabs(uu)>pp) {
        gto_int_2e_e_parm(&(ab[ip1]),&(ab[ip2]),pq,&th,&tt);
        /* coulombic integral */
        ee_int = gto_int_2e_e_core(&(ab[ip1]),&(ab[ip2]),pq,uu,th,tt,zr);
        gj_tot += (ip1==ip2 ? ee_int : 2.0*ee_int);
        /* exchange integral */
        gto_pair_swap(&ba,&(ab[ip2]));
        ee_int = gto_int_2e_e_core(&(ab[ip1]),&ba,pq,uu,th,tt,zr);
        gk_tot += (ip1==ip2 ? ee_int : 2.0*ee_int);
        }
      }
  /* clean memory */
  gto_pair_free(ab);
  /* finish */
  (*gj) = gj_tot;
  (*gk) = gk_tot;
  }

/* -------------------------------------------------------------------------- */
