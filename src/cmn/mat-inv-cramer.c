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
#include "cmn/matrix.h"
#include "cmn/message.h"

/* -------------------------------------------------------------------------- */

/* calculate inverse matrix by Cramer rule (double-precision real)
 
   m - the inverted matrix (output)
   u - the square matrix (input)
   n - dimension of the matrix */
void mat_finv_cramer(double **m, double **u, unsigned n) {
  unsigned i,j,k,l,i1,i2;
  double d,**w;
  /* analytical solution for small matrices */
  if (n<4)
    mat_finv(m,u,n);
  /* Cramer rule */
  else {
    w = mat_falloc(n-1,n-1);
    /* determinant */
    d = mat_fdet_cramer(u,n);
    /* adjugate matrix */
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
        /* submatrix */
        i1 = 0;
        for (k=0; k<n; k++) {
          if (k!=i) {
            i2 = 0;
            for (l=0; l<n; l++) {
              if (l!=j) {
                w[i1][i2] = u[k][l];
                i2++;
                }
              }
            i1++;
            }
          }
        /* inversion matrix element */
        m[j][i] = ((i+j)%2 ? -1.0 : 1.0)*mat_fdet_cramer(w,n-1)/d;
        }
      }
    mat_ffree(w,n-1);
    }
  }

/* calculate inverse matrix by Cramer rule (double-precision complex)
 
   m - the inverted matrix (output)
   u - the square matrix (input)
   n - dimension of the matrix */
void mat_zinv_cramer(double complex **m, double complex **u, unsigned n) {
  double complex d;
  double complex **w;
  unsigned i,j,k,l,i1,i2;
  /* analytical solution for small matrices */
  if (n<4)
    mat_zinv(m,u,n);
  /* Cramer rule */
  else {
    w = mat_zalloc(n-1,n-1);
    /* determinant */
    d = mat_zdet_cramer(u,n);
    /* adjugate matrix */
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
        /* submatrix */
        i1 = 0;
        for (k=0; k<n; k++) {
          if (k!=i) {
            i2 = 0;
            for (l=0; l<n; l++) {
              if (l!=j) {
                w[i1][i2] = u[k][l];
                i2++;
                }
              }
            i1++;
            }
          }
        /* inversion matrix element */
        m[j][i] = ((i+j)%2 ? -1.0 : 1.0)*mat_zdet_cramer(w,n-1)/d;
        }
      }
    mat_zfree(w,n-1);
    }
  }

/* -------------------------------------------------------------------------- */
