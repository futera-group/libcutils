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

/* calculate matrix determinant by Cramer rule (double-precision real)

   m - the matrix
   n - dimension of the matrix */
double mat_fdet_cramer(double **m, unsigned n) {
  unsigned i,j,k,id;
  double **w,d = 0.0;
  /* analytical solution for small matrices */
  if (n<4)
    d = mat_fdet(m,n);
  /* Cramer rule */
  else {
    w = mat_falloc(n-1,n-1);
    for (i=0; i<n; i++) {
      id = 0; 
      for (j=0; j<n; j++) {
        if (j!=i) {
          for (k=1; k<n; k++)
            w[id][k-1] = m[j][k];
          id++;
          }
        }
      d += (i%2 ? -1.0 : 1.0)*m[i][0]*mat_fdet_cramer(w,n-1);
      } 
    mat_ffree(w,n-1);
    }
  return(d);
  }

/* calculate matrix determinant by Cramer rule (double-precision complex)

   m - the matrix
   n - dimension of the matrix */
double complex mat_zdet_cramer(double complex **m, unsigned n) {
  double complex **w,d = 0.0;
  unsigned i,j,k,id;
  /* analytical solution for small matrices */
  if (n<4)
    d = mat_zdet(m,n);
  /* Cramer rule */
  else {
    w = mat_zalloc(n-1,n-1);
    for (i=0; i<n; i++) {
      id = 0; 
      for (j=0; j<n; j++) {
        if (j!=i) {
          for (k=1; k<n; k++)
            w[id][k-1] = m[j][k];
          id++;
          }
        }
      d += (i%2 ? -1.0 : 1.0)*m[i][0]*mat_zdet_cramer(w,n-1);
      } 
    mat_zfree(w,n-1);
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
