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
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

#define ORTH_CHECK_PREC 1.0E-10

/* -------------------------------------------------------------------------- */

/* Lowdin orthogonalization of m n-dimensional vectors

   a - matrix containing the vectors (in columns)
   n - dimension of the vectors
   m - number of the vectors */
void math_orth_lowdin_norm(double **a, unsigned n, unsigned m) {
  double **s,**u,**t,**x,*e;
  unsigned i,j;
  /* normalization */
  for (i=0; i<m; i++)
    vec_funit(a[i],n);
  /* overlap matrix */
  s = mat_falloc(m,m);
  for (i=0; i<m; i++)
    for (j=i; j<m; j++)
      s[i][j] = s[j][i] = vec_fprod_scalar(a[i],a[j],n);
  /* Lowdin transformation matrix */
  e = vec_falloc(m);
  u = mat_falloc(m,m);
  mat_diag(s,u,e,m,0);
  t = mat_falloc(m,m);
  for (i=0; i<m; i++) {
    t[i][i] = 1.0/sqrt(e[i]);
    for (j=i+1; j<m; j++)
      t[i][j] = 0.0;
    }
  x = mat_falloc(m,m);
  mat_ftrans(u,m);
  mat_fmult(x,t,u,m,m,m);
  mat_ftrans(u,m);
  mat_fmult(t,u,x,m,m,m);
  /* orthogonalization */
  mat_ftrans(a,m);
  for (i=0; i<m; i++) {
    mat_fmult_lvec(e,t,a[i],m,m);
    vec_fcopy(a[i],e,m);
    }
  mat_ftrans(a,m);
  /* clean memory */
  mat_ffree(s,m);
  mat_ffree(u,m);
  mat_ffree(t,m);
  mat_ffree(x,m);
  vec_ffree(e);
  }

/* -------------------------------------------------------------------------- */

/* Gramm-Schmidth orthogonalization of m n-dimensional vectors

   a - matrix containing the vectors (in columns)
   n - dimension of the vectors
   m - number of the vectors */
void math_orth_gs(double **a, unsigned n, unsigned m) {
  unsigned i,j,k;
  double s,*p,*u;
  if (m<2)
    return;
  /* allocate auxiliary vectors */
  p = vec_falloc(m);
  u = vec_falloc(m);
  /* calculate square of first vector */
  u[0] = 0.0;
  for (i=0; i<n; i++)
    u[0] += (a[i][0]*a[i][0]);
  /* keep first vector and orthogonalize the others */
  for (i=1; i<m; i++) {
    /* calculate projections */
    for (j=0; j<i; j++) {
      s = 0.0;
      for (k=0; k<n; k++)
        s += (a[k][j]*a[k][i]);
      p[j] = s/u[j];
      }
    /* orthogonalize vector */
    for (j=0; j<i; j++)
      for (k=0; k<n; k++)
        a[k][i] -= (p[j]*a[k][j]);
    /* calculate square of the vector */
    u[i] = 0.0;
    for (j=0; j<n; j++)
      u[i] += (a[j][i]*a[j][i]);
    }
  /* free allocated memory */
  vec_ffree(p);
  vec_ffree(u);
  }

/* Gramm-Schmidth orthonormalization of m n-dimensional vectors

   a - matrix containing the vectors (in columns)
   n - dimension of the vectors
   m - number of the vectors */
void math_orth_gs_norm(double **a, unsigned n, unsigned m) {
  unsigned i,j,k;
  double s,*p;
  if (m<2)
    return;
  /* allocate auxiliary vectors */
  p = vec_falloc(m);
  /* normalize first vector */
  s = 0.0;
  for (i=0; i<n; i++)
    s += (a[i][0]*a[i][0]);
  s = sqrt(s);
  for (i=0; i<n; i++)
    a[i][0] /= s;
  /* keep first vector and orthonormalize the others */
  for (i=1; i<m; i++) {
    /* calculate projections */
    for (j=0; j<i; j++) {
      p[j] = 0.0;
      for (k=0; k<n; k++)
        p[j] += (a[k][j]*a[k][i]);
      }
    /* orthogonalize vector */
    for (j=0; j<i; j++)
      for (k=0; k<n; k++)
        a[k][i] -= (p[j]*a[k][j]);
    /* normalize the vector */
    s = 0.0;
    for (j=0; j<n; j++)
      s += (a[j][i]*a[j][i]);
    s = sqrt(s);
    for (j=0; j<n; j++)
      a[j][i] /= s;
    }
  /* free allocated memory */
  vec_ffree(p);
  }

/* check orthogonality of given vectors

   a - matrix containing the vectors (in columns)
   n - dimension of the vectors
   m - number of the vectors */
short math_orth_check(double **a, unsigned n, unsigned m) {
  unsigned i,j,k;
  double s;
  for (i=0; i<m; i++)
    for (j=i; j<m; j++) {
      s = 0.0;
      for (k=0; k<n; k++)
        s += a[k][i]*a[k][j];
      if (i!=j && fabs(s)>ORTH_CHECK_PREC)
        return(0);
      }
  return(1);
  }

/* check orthonormality of given vectors

   a - matrix containing the vectors (in columns)
   n - dimension of the vectors
   m - number of the vectors */
short math_orth_check_norm(double **a, unsigned n, unsigned m) {
  unsigned i,j,k;
  double s;
  for (i=0; i<m; i++)
    for (j=i; j<m; j++) {
      s = 0.0;
      for (k=0; k<n; k++)
        s += a[k][i]*a[k][j];
      if (i==j && fabs(s-1.0)>ORTH_CHECK_PREC)
        return(0);
      if (i!=j && fabs(s)>ORTH_CHECK_PREC)
        return(0);
      }
  return(1);
  }

/* print table of vector products and return whether the vectors are orthogonal

   a - matrix containing the vectors (in columns)
   n - dimension of the vectors
   m - number of the vectors */
short math_orth_print(double **a, unsigned n, unsigned m) {
  short ret=1;
  unsigned i,j,k,mm;
  double s;
  mm = str_unum_len(m);
  for (i=0; i<m; i++)
    for (j=i; j<m; j++) {
      s = 0.0;
      for (k=0; k<n; k++)
        s += a[k][i]*a[k][j];
      printf("   (%*d,%*d) = %15.6E",mm,i+1,mm,j+1,s);
      if (i!=j && fabs(s)>ORTH_CHECK_PREC) {
        printf("  F");
        ret = 0;
        }
      printf("\n");
      }
  return(ret);
  }

/* print table of vector products and return if the vectors are orthonormal

   a - matrix containing the vectors (in columns)
   n - dimension of the vectors
   m - number of the vectors */
short math_orth_print_norm(double **a, unsigned n, unsigned m) {
  short ret=1;
  unsigned i,j,k,mm;
  double s;
  mm = str_unum_len(m);
  for (i=0; i<m; i++)
    for (j=i; j<m; j++) {
      s = 0.0;
      for (k=0; k<n; k++)
        s += a[k][i]*a[k][j];
      printf("   (%*d,%*d) = %15.6E",mm,i+1,mm,j+1,s);
      if (i==j && fabs(s-1.0)>ORTH_CHECK_PREC) {
        printf("  F");
        ret = 0;
        }
      if (i!=j && fabs(s)>ORTH_CHECK_PREC) {
        printf("  F");
        ret = 0;
        }
      printf("\n");
      }
  return(ret);
  }

/* -------------------------------------------------------------------------- */
