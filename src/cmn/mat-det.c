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

/* calculate matrix determinant up to size 3 (double-precision real)

   m - the matrix
   n - dimension of the matrix */
double mat_fdet(double **m, unsigned n) {
  double d = 0.0;
  /* dimension */
  switch (n) {
    /* none */
    case 0:
      d = 0.0;
      break;
    /* skalar */
    case 1: 
      d = m[0][0];
      break;
    /* 2x2 matrix */
    case 2:
      d = m[0][0]*m[1][1]-m[0][1]*m[1][0];
      break;
    /* 3x3 matrix */
    case 3:
      d  = m[0][0]*m[1][1]*m[2][2];
      d += m[0][2]*m[1][0]*m[2][1];
      d += m[0][1]*m[1][2]*m[2][0];
      d -= m[0][2]*m[1][1]*m[2][0];
      d -= m[0][1]*m[1][0]*m[2][2];
      d -= m[0][0]*m[1][2]*m[2][1];
      break;
    /* higher-order matrix */
    default:
      msg_error_f("determinant of matrix larger than 3x3",1);
      break;
    }
  return(d);
  }

/* calculate matrix determinant up to size 3 (double-precision complex)

   m - the matrix
   n - dimension of the matrix */
complex double mat_zdet(complex double **m, unsigned n) {
  double complex d = 0.0;
  /* dimension */
  switch (n) {
    /* none */
    case 0:
      d = 0.0;
      break;
    /* skalar */
    case 1:
      d = m[0][0];
      break;
    /* 2x2 matrix */
    case 2:
      d = m[0][0]*m[1][1]-m[0][1]*m[1][0];
      break;
    /* 3x3 matrix */
    case 3:
      d  = m[0][0]*m[1][1]*m[2][2];
      d += m[0][2]*m[1][0]*m[2][1];
      d += m[0][1]*m[1][2]*m[2][0];
      d -= m[0][2]*m[1][1]*m[2][0];
      d -= m[0][1]*m[1][0]*m[2][2];
      d -= m[0][0]*m[1][2]*m[2][1];
      break;
    /* higher-order matrix */
    default:
      msg_error_f("determinant of matrix larger than 3x3",1);
      break;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
