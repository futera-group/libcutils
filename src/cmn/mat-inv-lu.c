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

#include <complex.h>

/* -------------------------------------------------------------------------- */

/* calculate inverse matrix by LU decomposition (double-precision real)
 
   m - the inverted matrix (output)
   u - the LU-decomposed matrix (input)
   p - permutation index array from LU decomposition
   n - dimension of the matrix */
void mat_finv_lu(double **m, double **u, unsigned *p, unsigned n) {
  int i,j,k;
  for (j=0; j<n; j++) {
    for (i=0; i<n; i++) {
      m[i][j] = (p[i]==j ? 1.0 : 0.0);
      for (k=0; k<i; k++)
        m[i][j] -= u[i][k] * m[k][j];
      }
    for (i=n-1; i>=0; i--) {
      for (k=i+1; k<n; k++)
        m[i][j] -= u[i][k] * m[k][j];
      m[i][j] /= u[i][i];
      }
    }
  }

/* calculate inverse matrix by LU decomposition (double-precision complex)
 
   m - the inverted matrix (output)
   u - the LU-decomposed matrix (input)
   p - permutation index array from LU decomposition
   n - dimension of the matrix */
void mat_zinv_lu(complex double **m, complex double **u,
  unsigned *p, unsigned n) {
  int i,j,k;
  for (j=0; j<n; j++) {
    for (i=0; i<n; i++) {
      m[i][j] = (p[i]==j ? 1.0 : 0.0);
      for (k=0; k<i; k++)
        m[i][j] -= u[i][k] * m[k][j];
      }
    for (i=n-1; i>=0; i--) {
      for (k=i+1; k<n; k++)
        m[i][j] -= u[i][k] * m[k][j];
      m[i][j] /= u[i][i];
      }
    }
  }

/* -------------------------------------------------------------------------- */
