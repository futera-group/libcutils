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
#include "qmc/basis.h"
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* auxiliary function for printing value of 2e integrals
 
   s1,s2,s3,s4 - basis set shells
   c1,c2,c3,c4 - basis set shell coordinates
   f1,f2,f3,f4 - basis set function IDs
   up          - criterion for neglecting pair contributions
   uu          - criterion for neglecting pair-pair contributions
   zr          - criterion for neglecting the contribution
   f           - pointer to open output stream */
void gto_int_2e_print(
  struct basis_shell *s1, double *c1, unsigned f1,
  struct basis_shell *s2, double *c2, unsigned f2,
  struct basis_shell *s3, double *c3, unsigned f3,
  struct basis_shell *s4, double *c4, unsigned f4,
  double up, double uu, double zr, FILE *f) {
  static double v[2],r12[3],r34[3],r2_12,r2_34;
  static unsigned i,id[2][4],n = 0;
  static long unsigned n_tot = 0;
  if (s1 && s2 && s3 && s4) {
    /* intershell distances */
    r2_12 = r2_34 = 0.0;
    for (i=0; i<3; i++) {
      r12[i] = c1[i]-c2[i];
      r2_12 += r12[i]*r12[i];
      r34[i] = c3[i]-c4[i];
      r2_34 += r34[i]*r34[i];
      }
    /* (ab|cd) integral */
    v[n] = gto_int_2e_e_bf(s1,c1,f1,s2,c2,f2,s3,c3,f3,
      s4,c4,f4,r12,r2_12,r34,r2_34,up,uu,zr);
    /* save integral IDs */
    id[n][0] = s1->bf+f1;
    id[n][1] = s2->bf+f2;
    id[n][2] = s3->bf+f3;
    id[n][3] = s4->bf+f4;
    n++;
    /* print out */
    if (n>=2) {
      printf("%5ld: (%2d%2d|%2d%2d) = %14.6e %5ld: (%2d%2d|%2d%2d) = %14.6e\n",
        (n_tot+1)%10000,id[0][0]+1,id[0][1]+1,id[0][2]+1,id[0][3]+1,v[0],
        (n_tot+2)%10000,id[1][0]+1,id[1][1]+1,id[1][2]+1,id[1][3]+1,v[1]);
      n_tot += 2;
      n = 0;
      }
    }
  /* print out last integral */
  else if (n)
    printf("%5ld: (%2d%2d|%2d%2d) = %14.6e\n",
      (n_tot+1)%10000,id[0][0]+1,id[0][1]+1,id[0][2]+1,id[0][3]+1,v[0]);
  }

/* print out array with values of all symmetrically unique 2e AO integrals
 
   c1,c2   - basis set centers
   s1,s2   - basis set shells
   if1,if2 - basis shell function IDs
   id1,id2 - basis function IDs
   up      - criterion for neglecting pair contributions
   uu      - criterion for neglecting pair-pair contributions
   zr      - criterion for neglecting the contribution
   f       - pointer to open output stream */
void gto_int_2e_print_all_34(struct basis *b,
  struct basis_center *c1, struct basis_shell *s1, unsigned if1, unsigned id1,
  struct basis_center *c2, struct basis_shell *s2, unsigned if2, unsigned id2,
  double up, double uu, double zr, FILE *f) {
  unsigned ic3,ic4,is3,is4,if3,if4,id3,id4;
  struct basis_center *c3,*c4;
  struct basis_shell *s3,*s4;
  /* loop over third center */
  for (ic3=0; ic3<b->n_centers; ic3++) {
    c3 = b->center+ic3;
    for (is3=0; is3<c3->n_shells; is3++) {
      s3 = c3->shell+is3;
      for (if3=0; if3<s3->n_bfce; if3++) {
        id3 = s3->bf+if3;
        if (id1<id2 && id2<id3) {  
          gto_int_2e_print( /* (11|23) */
            s1,c1->coord,if1,s1,c1->coord,if1,
            s2,c2->coord,if2,s3,c3->coord,if3,up,uu,zr,f);
          gto_int_2e_print( /* (12|13) */
            s1,c1->coord,if1,s2,c2->coord,if2,
            s1,c1->coord,if1,s3,c3->coord,if3,up,uu,zr,f);
          }
        /* loop over fourth center */
        for (ic4=0; ic4<b->n_centers; ic4++) {
          c4 = b->center+ic4;
          for (is4=0; is4<c4->n_shells; is4++) {
            s4 = c4->shell+is4;
            for (if4=0; if4<s4->n_bfce; if4++) {
              id4 = s4->bf+if4;
              if (id1<id2 && id2<id3 && id3<id4) {
                gto_int_2e_print( /* (12|34) */
                  s1,c1->coord,if1,s2,c2->coord,if2,
                  s3,c3->coord,if3,s4,c4->coord,if4,up,uu,zr,f);
                gto_int_2e_print( /* (13|24) */
                  s1,c1->coord,if1,s3,c3->coord,if3,
                  s2,c2->coord,if2,s4,c4->coord,if4,up,uu,zr,f);
                gto_int_2e_print( /* (14|23) */
                  s1,c1->coord,if1,s4,c4->coord,if4,
                  s2,c2->coord,if2,s3,c3->coord,if3,up,uu,zr,f);
                }
              }
            }
          }
        }
      }
    }
  }

/* print out array with values of all symmetrically unique 2e AO integrals
 
   b  - pointer to basis set data
   up - criterion for neglecting pair contributions
   uu - criterion for neglecting pair-pair contributions
   zr - criterion for neglecting the contribution
   f  - pointer to open output stream */
void gto_int_2e_print_all(struct basis *b, double up, double uu,
  double zr, FILE *f) {
  unsigned ic1,ic2,is1,is2,if1,if2,id1,id2;
  struct basis_center *c1,*c2;
  struct basis_shell *s1,*s2;
  /* loop over first center */
  for (ic1=0; ic1<b->n_centers; ic1++) {
    c1 = b->center+ic1;
    for (is1=0; is1<c1->n_shells; is1++) {
      s1 = c1->shell+is1;
      for (if1=0; if1<s1->n_bfce; if1++) {
        id1 = s1->bf+if1;
        /* one-center integral */
        gto_int_2e_print( /* (11|11) */
          s1,c1->coord,if1,s1,c1->coord,if1,
          s1,c1->coord,if1,s1,c1->coord,if1,up,uu,zr,f);
        /* loop over second center */
        for (ic2=0; ic2<b->n_centers; ic2++) {
          c2 = b->center+ic2;
          for (is2=0; is2<c2->n_shells; is2++) {
            s2 = c2->shell+is2;
            for (if2=0; if2<s2->n_bfce; if2++) {
              id2 = s2->bf+if2;
              if (id1!=id2)
                gto_int_2e_print( /* (11|12) */
                  s1,c1->coord,if1,s1,c1->coord,if1,
                  s1,c1->coord,if1,s2,c2->coord,if2,up,uu,zr,f);
              if (id1<id2) { 
                gto_int_2e_print( /* (11|22) */
                  s1,c1->coord,if1,s1,c1->coord,if1,
                  s2,c2->coord,if2,s2,c2->coord,if2,up,uu,zr,f);
                gto_int_2e_print( /* (12|12) */
                  s1,c1->coord,if1,s2,c2->coord,if2,
                  s1,c1->coord,if1,s2,c2->coord,if2,up,uu,zr,f);
                }
              /* loop over third and fourth center */
              gto_int_2e_print_all_34(b,c1,s1,if1,id1,c2,s2,if2,id2,up,uu,zr,f);
              }
            }
          }
        }
      }
    }
  gto_int_2e_print( /* (**|**) */
    NULL,NULL,0,NULL,NULL,0,
    NULL,NULL,0,NULL,NULL,0,up,uu,zr,f);
  }

/* -------------------------------------------------------------------------- */
