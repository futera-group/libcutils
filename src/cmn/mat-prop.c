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
#include "cmn/matrix.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* check if two double real matrices are the same

   m - double precision real square matrix
   n - dimension of the matrix */
short mat_fis_same(double **a, double **b, unsigned n, unsigned m) {
  unsigned i,j;
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      if (fabs(a[i][j]-b[i][j])>MAT_CHECK_PREC)
        return(0);
  return(1);
  }

/* check if two double real matrices are the same and print differences

   m - double precision real square matrix
   n - dimension of the matrix */
short mat_fis_same_print(double **a, double **b, unsigned n, unsigned m) {
  short ret = 1;
  unsigned i,j,nn;
  double diff;
  nn = str_unum_len(n);
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      diff = fabs(a[i][j]-b[i][j]);
      printf("   (%*d,%*d): %15.6E",nn,i+1,nn,j+1,diff);
      if (diff>MAT_CHECK_PREC) {
        printf("  F");
        ret=0;
        }
      printf("\n");
      }
  return(ret);
  }

/* check if double real matrix is symmetric

   m - double precision real square matrix
   n - dimension of the matrix */
short mat_fis_sym(double **m, unsigned n) {
  unsigned i,j;
  for (i=0; i<n; i++)
    for (j=i; j<n; j++)
      if (fabs(m[i][j]-m[j][i])>MAT_CHECK_PREC)
        return(0);
  return(1);
  }

/* check if double real matrix is symmetric and print off-diag differences

   m - double precision real square matrix
   n - dimension of the matrix */
short mat_fis_sym_print(double **m, unsigned n) {
  short ret=1;
  unsigned i,j,nn;
  double diff;
  nn = str_unum_len(n);
  for (i=0; i<n; i++)
    for (j=i; j<n; j++) {
      diff = fabs(m[i][j]-m[j][i]);
      printf("   (%*d,%*d): %15.6E",nn,i+1,nn,j+1,diff);
      if (diff>MAT_CHECK_PREC) {
        printf("  F");
        ret=0;
        }
      printf("\n");
      }
  return(ret);
  }

/* check if double-precision real matrix is orthogonal
 
   m     - the matrix
   n1,n2 - dimensions of the matrix */
short mat_fis_orth(double **m, unsigned n1, unsigned n2) {
  unsigned i,j,k;
  double s;
  /* check row products */
  for (i=0; i<n1; i++)
    for (j=i; j<n1; j++) {
      s = 0.0;
      for (k=0; k<n2; k++)
        s += m[i][k]*m[j][k];
      if (i==j && (fabs(s-1.0)>MAT_CHECK_PREC))
        return(0);
      if (i!=j && (fabs(s)>MAT_CHECK_PREC))
        return(0);
      }
  return(1);
  }

/* check if double real matrix is orthogonal

   m - double precision real square matrix
   n - dimension of the matrix */
short mat_fis_orth_Sc(double **m, unsigned n) {
  unsigned i,j,k;
  double sc,sr;
  for (i=0; i<n; i++)
    for (j=i; j<n; j++) {
      sc = sr = 0.0;
      for (k=0; k<n; k++) {
        sc += m[i][k]*m[j][k];
        sr += m[k][i]*m[k][j];
        }
      if (i==j && (fabs(sc-1.0)>MAT_CHECK_PREC || fabs(sr-1.0)>MAT_CHECK_PREC))
        return(0);
      if (i!=j && (fabs(sc)>MAT_CHECK_PREC || fabs(sr)>MAT_CHECK_PREC))
        return(0);
      }
  return(1);
  }

/* check if double real matrix is orthogonal and print all vector products

   m - double precision real square matrix
   n - dimension of the matrix */
short mat_fis_orth_Sc_print(double **m, unsigned n) {
  short ret = 1;
  unsigned i,j,k,nn;
  double sc,sr;
  nn = str_unum_len(n);
  for (i=0; i<n; i++)
    for (j=i; j<n; j++) {
      sc = sr = 0.0;
      for (k=0; k<n; k++) {
        sr += m[i][k]*m[j][k];
        sc += m[k][i]*m[k][j];
        }
      printf("   (%*d,%*d): Col = %15.6E, Row = %15.6E",nn,i+1,nn,j+1,sc,sr);
      if (i==j) {
        if (fabs(sc-1.0)>MAT_CHECK_PREC)
          printf("  F(c)");
        if (fabs(sr-1.0)>MAT_CHECK_PREC)
          printf("  F(r)");
        ret = 0;
        }
      else {
        if (fabs(sc)>MAT_CHECK_PREC)
          printf("  F(c)");
        if (fabs(sr)>MAT_CHECK_PREC)
          printf("  F(r)");
        ret=0;
        }
      printf("\n");
      }
  return(ret);
  }

/* check if double real matrix is unit matrix

   m - double precision real square matrix
   n - dimension of the matrix
   e - largest discrepancy (optional) */
short mat_fis_unit(double **m, unsigned n, double *e) {
  unsigned i,j;
  double err = 0.0;
  for (i=0; i<n; i++) {
    if (fabs(m[i][i]-1.0)>err)
      err = fabs(m[i][i]-1.0);
    for (j=i+1; j<n; j++) {
      if (fabs(m[i][j])>err) 
        err = fabs(m[i][j]);
      if (fabs(m[j][i])>err) 
        err = fabs(m[j][i]);
      }
    }
  if (*e)
    (*e) = err;
  return(err<=MAT_CHECK_PREC);
  }

/* check if double real matrix is diagonal

   m - double precision real square matrix
   n - dimension of the matrix */
short mat_fis_diag(double **m, unsigned n) {
  unsigned i,j;
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++)
      if (fabs(m[i][j])>MAT_CHECK_PREC)
        return(0);
  return(1);
  }

/* -------------------------------------------------------------------------- */
