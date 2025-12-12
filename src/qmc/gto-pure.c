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

#include <cmn/math.h>
#include <cmn/message.h>
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* transform D5 GTO to linear combination of Cartesian GTO
 
   s     - sub-type of spherical GTO (p-x,p-y,p-z,...)
   lc_id - array for Cartesian GTO IDs (3 per GTO = l,m,n)
   lc_mf - array for multiplication factors in the linear combination
   lc_n  - number of terms in the linear combination */
void gto_pure_to_cart_d5(short s, unsigned *lc_id,
  double *lc_mf, unsigned *lc_n) {
  switch (s) {
    case 0: /* T(2,-2) = D(xy) */
      lc_id[0] = 1; 
      lc_id[1] = 1;
      lc_id[2] = 0;
      lc_mf[0] = 1.0;
      (*lc_n) = 1;
      break;
    case 1: /* T(2,-1) = D(yz) */
      lc_id[0] = 0;
      lc_id[1] = 1;
      lc_id[2] = 1;
      lc_mf[0] = 1.0;
      (*lc_n) = 1;
      break;
    case 2: /* T(2, 0) = D(z2) */
      lc_id[0] = 0;
      lc_id[1] = 0;
      lc_id[2] = 2;
      lc_id[3] = 2;
      lc_id[4] = 0; 
      lc_id[5] = 0;
      lc_id[6] = 0;
      lc_id[7] = 2;
      lc_id[8] = 0;
      lc_mf[0] = 1.0;
      lc_mf[1] = -0.5;
      lc_mf[2] = -0.5;
      (*lc_n) = 3;
      break;
    case 3: /* T(2,+1) = D(xz) */
      lc_id[0] = 1;
      lc_id[1] = 0;
      lc_id[2] = 1;
      lc_mf[0] = 1.0;
      (*lc_n) = 1;
      break;
    case 4: /* T(2,+2) = D(x2-y2) */
      lc_id[0] = 2;
      lc_id[1] = 0;
      lc_id[2] = 0;
      lc_id[3] = 0;
      lc_id[4] = 2;
      lc_id[5] = 0;
      lc_mf[0] = MATH_2r3f2;
      lc_mf[1] = -MATH_2r3f2;
      (*lc_n) = 2;
      break;
    }
  }

/* transform F7 GTO to linear combination of Cartesian GTO
 
   s     - sub-type of spherical GTO (p-x,p-y,p-z,...)
   lc_id - array for Cartesian GTO IDs (3 per GTO = l,m,n)
   lc_mf - array for multiplication factors in the linear combination
   lc_n  - number of terms in the linear combination */
void gto_pure_to_cart_f7(short s, unsigned *lc_id,
  double *lc_mf, unsigned *lc_n) {
  switch (s) {
    case 0: /* T(3,-3) = F(x2y-y3) */
      lc_id[0] = 2;
      lc_id[1] = 1;
      lc_id[2] = 0;
      lc_id[3] = 0;
      lc_id[4] = 3;
      lc_id[5] = 0;
      lc_mf[0] = MATH_3f2m2r2;
      lc_mf[1] = -MATH_2r5f2m2r2;
      (*lc_n) = 2;
      break;
      break;
    case 1: /* T(3,-2) = F(xyz) */
      lc_id[0] = 1;
      lc_id[1] = 1;
      lc_id[2] = 1;
      lc_mf[0] = 1.0;
      (*lc_n) = 1;
      break;
    case 2: /* T(3,-1) = F(z2y) */
      lc_id[0] = 0;
      lc_id[1] = 1;
      lc_id[2] = 2;
      lc_id[3] = 2;
      lc_id[4] = 1;
      lc_id[5] = 0;
      lc_id[6] = 0;
      lc_id[7] = 3;
      lc_id[8] = 0;
      lc_mf[0] = MATH_2r6f2r5;
      lc_mf[1] = -MATH_2r3f2m2r10;
      lc_mf[2] = -MATH_2r3f2m2r2;
      (*lc_n) = 3;
      break;
    case 3: /* T(3, 0) = F(z3) */
      lc_id[0] = 0;
      lc_id[1] = 0;
      lc_id[2] = 3;
      lc_id[3] = 0;
      lc_id[4] = 2;
      lc_id[5] = 1;
      lc_id[6] = 2;
      lc_id[7] = 0;
      lc_id[8] = 1;
      lc_mf[0] = 1.0;
      lc_mf[1] = -MATH_2r3f2m2r5;
      lc_mf[2] = -MATH_2r3f2m2r5;
      (*lc_n) = 3;
      break;
    case 4: /* T(3,+1) = F(z2x) */
      lc_id[0] = 1;
      lc_id[1] = 0;
      lc_id[2] = 2;
      lc_id[3] = 3;
      lc_id[4] = 0;
      lc_id[5] = 0;
      lc_id[6] = 1;
      lc_id[7] = 2;
      lc_id[8] = 0;
      lc_mf[0] = MATH_2r6f2r5;
      lc_mf[1] = -MATH_2r3f2m2r2;
      lc_mf[2] = -MATH_2r3f2m2r10;
      (*lc_n) = 3;
      break;
    case 5: /* T(3,+2) = F(x2z-y2z) */
      lc_id[0] = 2;
      lc_id[1] = 0;
      lc_id[2] = 1;
      lc_id[3] = 0;
      lc_id[4] = 2;
      lc_id[5] = 1;
      lc_mf[0] = MATH_2r3f2;
      lc_mf[1] = -MATH_2r3f2;
      (*lc_n) = 2;
      break;
    case 6: /* T(3,+3) = F(x3-y2x) */
      lc_id[0] = 3;
      lc_id[1] = 0;
      lc_id[2] = 0;
      lc_id[3] = 1;
      lc_id[4] = 2;
      lc_id[5] = 0;
      lc_mf[0] = MATH_2r5f2m2r2;
      lc_mf[1] = -MATH_3f2m2r2;
      (*lc_n) = 2;
      break;
    }
  }

/* transform G9 GTO to linear combination of Cartesian GTO
 
   s     - sub-type of spherical GTO (p-x,p-y,p-z,...)
   lc_id - array for Cartesian GTO IDs (3 per GTO = l,m,n)
   lc_mf - array for multiplication factors in the linear combination
   lc_n  - number of terms in the linear combination */
void gto_pure_to_cart_g9(short s, unsigned *lc_id,
  double *lc_mf, unsigned *lc_n) {
  switch (s) {
    case 0: /* T(4,-4) = G(x3y-xy3) */
      lc_id[0] = 3;
      lc_id[1] = 1;
      lc_id[2] = 0;
      lc_id[3] = 1;
      lc_id[4] = 3;
      lc_id[5] = 0;
      lc_mf[0] = MATH_2r5f8;
      lc_mf[1] = -MATH_2r5f8;
      (*lc_n) = 2;
      break;
    case 1: /* T(4,-3) = G(x2yz-y3z) */
      lc_id[0] = 2;
      lc_id[1] = 1;
      lc_id[2] = 1;
      lc_id[3] = 0;
      lc_id[4] = 3;
      lc_id[5] = 1;
      lc_mf[0] = MATH_3f2m2r2;
      lc_mf[1] = -MATH_2r5f2m2r2;
      (*lc_n) = 2;
      break;
    case 2: /* T(4,-2) = G(xyz2) */
      lc_id[0] = 1;
      lc_id[1] = 1;
      lc_id[2] = 2;
      lc_id[3] = 3;
      lc_id[4] = 1;
      lc_id[5] = 0;
      lc_id[6] = 1;
      lc_id[7] = 3;
      lc_id[8] = 0;
      lc_mf[0] = MATH_3f2r7;
      lc_mf[1] = -MATH_2r5f2m2r7;
      lc_mf[2] = -MATH_2r5f2m2r7;
      (*lc_n) = 3;
      break;
    case 3: /* T(4,-1) = G(yz3) */
      lc_id[0] = 0;
      lc_id[1] = 1;
      lc_id[2] = 3;
      lc_id[3] = 0;
      lc_id[4] = 3;
      lc_id[5] = 1;
      lc_id[6] = 2;
      lc_id[7] = 1;
      lc_id[8] = 1;
      lc_mf[0] = MATH_2r10f2r7;
      lc_mf[1] = -MATH_3m2r5f2m2r14;
      lc_mf[2] = -MATH_3f2m2r14;
      (*lc_n) = 3;
      break;
    case 4: /* T(4, 0) = G(z4) */
      lc_id[0] = 0;
      lc_id[1] = 0;
      lc_id[2] = 4;
      lc_id[3] = 4;
      lc_id[4] = 0;
      lc_id[5] = 0;
      lc_id[6] = 0;
      lc_id[7] = 4;
      lc_id[8] = 0;
      lc_id[9] = 2;
      lc_id[10] = 0;
      lc_id[11] = 2;
      lc_id[12] = 0;
      lc_id[13] = 2;
      lc_id[14] = 2;
      lc_id[15] = 2;
      lc_id[16] = 2;
      lc_id[17] = 0;
      lc_mf[0] = 1.0;
      lc_mf[1] = 0.375;
      lc_mf[2] = 0.375;
      lc_mf[3] = -MATH_3m2r3f2r35;
      lc_mf[4] = -MATH_3m2r3f2r35;
      lc_mf[5] = MATH_3m2r3f4m2r35;
      (*lc_n) = 6;
      break;
    case 5: /* T(4,+1) = G(xz3) */
      lc_id[0] = 1;
      lc_id[1] = 0;
      lc_id[2] = 3;
      lc_id[3] = 3;
      lc_id[4] = 0;
      lc_id[5] = 1;
      lc_id[6] = 1;
      lc_id[7] = 2;
      lc_id[8] = 1;
      lc_mf[0] = MATH_2r10f2r7;
      lc_mf[1] = -MATH_3m2r5f2m2r14;
      lc_mf[2] = -MATH_3f2m2r14;
      (*lc_n) = 3;
      break;
    case 6: /* T(4,+2) = G(x2z2-y2z2) */
      lc_id[0] = 2;
      lc_id[1] = 0;
      lc_id[2] = 2;
      lc_id[3] = 0;
      lc_id[4] = 2;
      lc_id[5] = 2;
      lc_id[6] = 4;
      lc_id[7] = 0;
      lc_id[8] = 0;
      lc_id[9] = 0;
      lc_id[10] = 4;
      lc_id[11] = 0;
      lc_mf[0] = MATH_3m2r3f2m2r7;
      lc_mf[1] = -MATH_3m2r3f2m2r7;
      lc_mf[2] = -MATH_2r5f4;
      lc_mf[3] = MATH_2r5f4;
      (*lc_n) = 4;
      break;
    case 7: /* T(4,+3) = G(x3z-xy2z) */
      lc_id[0] = 3;
      lc_id[1] = 0; 
      lc_id[2] = 1;
      lc_id[3] = 1;
      lc_id[4] = 2;
      lc_id[5] = 1;
      lc_mf[0] = MATH_2r5f2m2r2;
      lc_mf[1] = -MATH_3f2m2r2;
      (*lc_n) = 2;
      break;
    case 8: /* T(4,+4) = G(x4-x2y2+y4) */
      lc_id[0] = 4;
      lc_id[1] = 0;
      lc_id[2] = 0;
      lc_id[3] = 0;
      lc_id[4] = 4;
      lc_id[5] = 0;
      lc_id[6] = 2;
      lc_id[7] = 2;
      lc_id[8] = 0;
      lc_mf[0] = MATH_2r35f8;
      lc_mf[1] = MATH_2r35f8;
      lc_mf[2] = -MATH_3m2r3f4;
      (*lc_n) = 3;
      break;
    }
  }

/* transform H11 GTO to linear combination of Cartesian GTO
 
   s     - sub-type of spherical GTO (p-x,p-y,p-z,...)
   lc_id - array for Cartesian GTO IDs (3 per GTO = l,m,n)
   lc_mf - array for multiplication factors in the linear combination
   lc_n  - number of terms in the linear combination */
void gto_pure_to_cart_h11(short s, unsigned *lc_id,
  double *lc_mf, unsigned *lc_n) {
  switch (s) {
    case  0: /* T(5,-5) = H(x5+x4y-x2y2) */
      lc_id[0] = 0;
      lc_id[1] = 5;
      lc_id[2] = 0;
      lc_id[3] = 4;
      lc_id[4] = 1;
      lc_id[5] = 0;
      lc_id[6] = 2;
      lc_id[7] = 3;
      lc_id[8] = 0;
      lc_mf[0] = MATH_3m2r7f8m2r2;
      lc_mf[1] = MATH_5m2r7f8m2r2;
      lc_mf[2] = -MATH_5m2r3f4m2r2;
      (*lc_n) = 3;
      break;
    case  1: /* T(5,-4) = H(x3yz-xy3z) */
      lc_id[0] = 3;
      lc_id[1] = 1;
      lc_id[2] = 1;
      lc_id[3] = 1;
      lc_id[4] = 3;
      lc_id[5] = 1;
      lc_mf[0] = MATH_2r5f2;
      lc_mf[1] = -MATH_2r5f2;
      (*lc_n) = 2;
      break;
    case  2: /* T(5,-3) = H(y5-z4y) */
      lc_id[0] = 2;
      lc_id[1] = 1;
      lc_id[2] = 2;
      lc_id[3] = 0;
      lc_id[4] = 3;
      lc_id[5] = 2;
      lc_id[6] = 0;
      lc_id[7] = 5;
      lc_id[8] = 0;
      lc_id[9] = 4;
      lc_id[10] = 1;
      lc_id[11] = 0;
      lc_id[12] = 2;
      lc_id[13] = 3; 
      lc_id[14] = 0;
      lc_mf[0] = MATH_2r3f2r2;
      lc_mf[1] = -MATH_2r5f2r6;
      lc_mf[2] = MATH_2r35f8m2r2;
      lc_mf[3] = -MATH_2r35f8m2r2;
      lc_mf[4] = -MATH_2r5f4m2r6;
      (*lc_n) = 5;
      break;
    case  3: /* T(5,-2) = H(xyz3-x3yz-xy3z) */
      lc_id[0] = 1;
      lc_id[1] = 1;
      lc_id[2] = 3;
      lc_id[3] = 3;
      lc_id[4] = 1;
      lc_id[5] = 1;
      lc_id[6] = 1;
      lc_id[7] = 3;
      lc_id[8] = 1;
      lc_mf[0] = MATH_2r5f2r3;
      lc_mf[1] = -MATH_2r5f2m2r3;
      lc_mf[2] = -MATH_2r5f2m2r3;
      (*lc_n) = 3;
      break;
    case  4: /* T(5,-1) = H(y5) */
      lc_id[0] = 0;
      lc_id[1] = 1;
      lc_id[2] = 4;
      lc_id[3] = 0;
      lc_id[4] = 3;
      lc_id[5] = 2;
      lc_id[6] = 2;
      lc_id[7] = 1;
      lc_id[8] = 2;
      lc_id[9] = 0;
      lc_id[10] = 5;
      lc_id[11] = 0;
      lc_id[12] = 4;
      lc_id[13] = 1;
      lc_id[14] = 0;
      lc_id[15] = 2;
      lc_id[16] = 3;
      lc_id[17] = 0;
      lc_mf[0] = MATH_2r5f2r3;
      lc_mf[1] = -MATH_3m2r5f2m2r7;
      lc_mf[2] = -MATH_3f2m2r7;
      lc_mf[3] = MATH_2r15f8;
      lc_mf[4] = MATH_2r5f8m2r3;
      lc_mf[5] = MATH_2r5f4m2r7;
      (*lc_n) = 6;
      break;
    case  5: /* T(5, 0) = H(z5) */
      lc_id[0] = 0;
      lc_id[1] = 0;
      lc_id[2] = 5;
      lc_id[3] = 2;
      lc_id[4] = 0;
      lc_id[5] = 3;
      lc_id[6] = 0;
      lc_id[7] = 2;
      lc_id[8] = 3;
      lc_id[9] = 4;
      lc_id[10] = 0;
      lc_id[11] = 1;
      lc_id[12] = 0;
      lc_id[13] = 4;
      lc_id[14] = 1;
      lc_id[15] = 2; 
      lc_id[16] = 2;
      lc_id[17] = 1;
      lc_mf[0] = 1.0;
      lc_mf[1] = -MATH_5f2r21;
      lc_mf[2] = -MATH_5f2r21;
      lc_mf[3] = 0.625;
      lc_mf[4] = 0.625;
      lc_mf[5] = MATH_2r15f4m2r7;
      (*lc_n) = 6;
      break;
    case  6: /* T(5,+1) = H(x5) */
      lc_id[0] = 1;
      lc_id[1] = 0;
      lc_id[2] = 4;
      lc_id[3] = 3;
      lc_id[4] = 0;
      lc_id[5] = 2;
      lc_id[6] = 1;
      lc_id[7] = 2; 
      lc_id[8] = 2;
      lc_id[9] = 5;
      lc_id[10] = 0;
      lc_id[11] = 0;
      lc_id[12] = 1;
      lc_id[13] = 4;
      lc_id[14] = 0;
      lc_id[15] = 3;
      lc_id[16] = 2;
      lc_id[17] = 0;
      lc_mf[0] = MATH_2r5f2r3;
      lc_mf[1] = -MATH_3m2r5f2m2r7;
      lc_mf[2] = -MATH_3f2m2r7;
      lc_mf[3] = MATH_2r15f8;
      lc_mf[4] = MATH_2r5f8m2r3;
      lc_mf[5] = MATH_2r5f4m2r7;
      (*lc_n) = 6;
      break;
    case  7: /* T(5,+2) = H(y4z-x4z) */
      lc_id[0] = 2;
      lc_id[1] = 0;
      lc_id[2] = 3;
      lc_id[3] = 0;
      lc_id[4] = 2;
      lc_id[5] = 3;
      lc_id[6] = 4;
      lc_id[7] = 0;
      lc_id[8] = 1;
      lc_id[9] = 0;
      lc_id[10] = 4;
      lc_id[11] = 1;
      lc_mf[0] = MATH_2r5f2;
      lc_mf[1] = -MATH_2r5f2;
      lc_mf[2] = -MATH_2r35f4m2r3;
      lc_mf[3] = MATH_2r35f4m2r3;
      (*lc_n) = 4;
      break;
    case  8: /* T(5,+3) = H(xy4-x5) */
      lc_id[0] = 3;
      lc_id[1] = 0;
      lc_id[2] = 2;
      lc_id[3] = 1;
      lc_id[4] = 2;
      lc_id[5] = 2;
      lc_id[6] = 5;
      lc_id[7] = 0;
      lc_id[8] = 0;
      lc_id[9] = 1;
      lc_id[10] = 4;
      lc_id[11] = 0;
      lc_id[12] = 3;
      lc_id[13] = 2;
      lc_id[14] = 0;
      lc_mf[0] = MATH_2r5f2r6;
      lc_mf[1] = -MATH_2r3f2r2;
      lc_mf[2] = -MATH_2r35f8m2r2;
      lc_mf[3] = MATH_2r35f8m2r2;
      lc_mf[4] = MATH_2r5f4m2r6;
      (*lc_n) = 5;
      break;
    case  9: /* T(5,+4) = H(x4z+y4z-x2y2z) */
      lc_id[0] = 4;
      lc_id[1] = 0;
      lc_id[2] = 1;
      lc_id[3] = 0;
      lc_id[4] = 4;
      lc_id[5] = 1;
      lc_id[6] = 2;
      lc_id[7] = 2;
      lc_id[8] = 1;
      lc_mf[0] = MATH_2r35f8;
      lc_mf[1] = MATH_2r35f8;
      lc_mf[2] = -MATH_3m2r3f4;
      (*lc_n) = 3;
      break;
    case 10: /* T(5,+5) = H(x5+xy4-x2y3) */
      lc_id[0] = 5;
      lc_id[1] = 0;
      lc_id[2] = 0;
      lc_id[3] = 1;
      lc_id[4] = 4;
      lc_id[5] = 0;
      lc_id[6] = 3;
      lc_id[7] = 2;
      lc_id[8] = 0;
      lc_mf[0] = MATH_3m2r7f8m2r2;
      lc_mf[1] = MATH_5m2r7f8m2r2;
      lc_mf[2] = -MATH_5m2r3f4m2r2;
      (*lc_n) = 3;
      break;
    }
  }

/* transform spherical GTO to linear combination of Cartesian GTO
 
   t     - type of spherical GTO (s,p,d,...)
   s     - sub-type of spherical GTO (p-x,p-y,p-z,...)
   lc_id - array for Cartesian GTO IDs (3 per GTO = l,m,n)
   lc_mf - array for multiplication factors in the linear combination
   lc_n  - number of terms in the linear combination */
void gto_pure_to_cart(short t, short s,
  unsigned *lc_id, double *lc_mf, unsigned *lc_n) {
  (*lc_n) = 0;
  /* divide SP shell */
  if (t==BASIS_SHELL_SP) {
    t = (s>0 ? BASIS_SHELL_P : BASIS_SHELL_S);
    s = (s>0 ? s-1 : s);
    }
  /* set LC index array */
  switch (t) {
    /* Cartesian GTO */
    case BASIS_SHELL_S:
    case BASIS_SHELL_P:
    case BASIS_SHELL_Dc:
    case BASIS_SHELL_Fc:
    case BASIS_SHELL_Gc:
    case BASIS_SHELL_Hc:
      gto_ang_vec(t,s,lc_id);
      lc_mf[0] = 1.0;
      (*lc_n) = 1;
      break;
    /* Spherical GTO */
    case BASIS_SHELL_Dp:
      gto_pure_to_cart_d5(s,lc_id,lc_mf,lc_n);
      break;
    case BASIS_SHELL_Fp:
      gto_pure_to_cart_f7(s,lc_id,lc_mf,lc_n);
      break;
    case BASIS_SHELL_Gp:
      gto_pure_to_cart_g9(s,lc_id,lc_mf,lc_n);
      break;
    case BASIS_SHELL_Hp:
      gto_pure_to_cart_h11(s,lc_id,lc_mf,lc_n);
      break;
    default:
      msg_error_f("unsupported GTO type for pure->cart transformation:"
        " %d/%d",1,t,s);
    }
  }

/* -------------------------------------------------------------------------- */
