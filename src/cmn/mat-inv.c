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

#include "cmn/matrix.h"
#include "cmn/message.h"

/* -------------------------------------------------------------------------- */

/* calculate inverse matrix up to size 3 (double-precision real)
 
   m - the inverted matrix (output)
   u - the square matrix (input)
   n - dimension of the matrix */
void mat_finv(double **m, double **u, unsigned n) {
  double d;
  switch (n) {
    /* none */
    case 0:
      break;
    /* skalar */
    case 1:
      m[0][0] = 1.0/u[0][0];
      break;
    /* 2x2 matrix */
    case 2:
      d = mat_fdet_cramer(u,n);
      m[0][0] =  u[1][1]/d;
      m[0][1] = -u[1][0]/d;
      m[1][0] = -u[0][1]/d;
      m[1][1] =  u[0][0]/d;
      break;
    /* 3x3 matrix */
    case 3:
      d = mat_fdet_cramer(u,n);
      m[0][0] =  (u[1][1]*u[2][2]-u[1][2]*u[2][1])/d;
      m[0][1] = -(u[0][1]*u[2][2]-u[0][2]*u[2][1])/d;
      m[0][2] =  (u[0][1]*u[1][2]-u[0][2]*u[1][1])/d;
      m[1][0] = -(u[1][0]*u[2][2]-u[1][2]*u[2][0])/d;
      m[1][1] =  (u[0][0]*u[2][2]-u[0][2]*u[2][0])/d;
      m[1][2] = -(u[0][0]*u[1][2]-u[0][2]*u[1][0])/d;
      m[2][0] =  (u[1][0]*u[2][1]-u[1][1]*u[2][0])/d;
      m[2][1] = -(u[0][0]*u[2][1]-u[0][1]*u[2][0])/d;
      m[2][2] =  (u[0][0]*u[1][1]-u[0][1]*u[1][0])/d;
      break;
    /* higher-order matrix */
    default:
      msg_error_f("inversion of matrix larger than 3x3",1);
      break;
    }
  }

/* calculate inverse matrix up to size 3 (double-precision complex)
 
   m - the inverted matrix (output)
   u - the square matrix (input)
   n - dimension of the matrix */
void mat_zinv(complex double **m, complex double **u, unsigned n) {
  double complex d;
  /* dimension */
  switch (n) {
    /* none */
    case 0:
      break;
    /* skalar */
    case 1:
      m[0][0] = 1.0/u[0][0];
      break;
    /* 2x2 matrix */
    case 2:
      d = mat_zdet_cramer(u,n);
      m[0][0] =  u[1][1]/d;
      m[0][1] = -u[1][0]/d;
      m[1][0] = -u[0][1]/d;
      m[1][1] =  u[0][0]/d;
      break;
    /* 3x3 matrix */
    case 3:
      d = mat_zdet_cramer(u,n);
      m[0][0] =  (u[1][1]*u[2][2]-u[1][2]*u[2][1])/d;
      m[0][1] = -(u[0][1]*u[2][2]-u[0][2]*u[2][1])/d;
      m[0][2] =  (u[0][1]*u[1][2]-u[0][2]*u[1][1])/d;
      m[1][0] = -(u[1][0]*u[2][2]-u[1][2]*u[2][0])/d;
      m[1][1] =  (u[0][0]*u[2][2]-u[0][2]*u[2][0])/d;
      m[1][2] = -(u[0][0]*u[1][2]-u[0][2]*u[1][0])/d;
      m[2][0] =  (u[1][0]*u[2][1]-u[1][1]*u[2][0])/d;
      m[2][1] = -(u[0][0]*u[2][1]-u[0][1]*u[2][0])/d;
      m[2][2] =  (u[0][0]*u[1][1]-u[0][1]*u[1][0])/d;
      break;
    /* higher-order matrix */
    default:
      msg_error_f("inversion of matrix larger than 3x3",1);
      break;
    }
  }

/* -------------------------------------------------------------------------- */
