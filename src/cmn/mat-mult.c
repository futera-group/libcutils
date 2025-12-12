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

#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* multiplication of two matrices (integer version)

   m        - result matrix of dimension n1 x n3
   m1       - first matrix of dimension n1 x n2
   m2       - second matrix of dimension n2 x n3
   n1,n2,n3 - dimenstions of matrices */
void mat_imult(int **m, int **m1, int **m2,
  unsigned n1, unsigned n2, unsigned n3) {
  unsigned i,j,k;
  for (i=0; i<n1; i++)
    for (j=0; j<n3; j++) {
      m[i][j] = 0;
      for (k=0; k<n2; k++)
        m[i][j] += (m1[i][k]*m2[k][j]);
      }
  }

/* multiplication of matrix by vector from left (integer version)

   r     - resulting vector of dimension n2
   m     - matrix of dimension n1 x n2
   v     - vector of dimension n1
   n1,n2 - dimensions of matrix and vectors */
void mat_imult_lvec(int *r, int **m, int *v, unsigned n1,
  unsigned n2) {
  unsigned i,j;
  for (i=0; i<n2; i++) {
    r[i] = 0;
    for (j=0; j<n1; j++)
      r[i] += (v[j]*m[j][i]);
    }
  }

/* multiplication of matrix by vector from right (integer version)

   r     - resulting vector of dimension n1
   m     - matrix of dimension n1 x n2
   v     - vector of dimension n2
   n1,n2 - dimensions of matrix and vectors */
void mat_imult_rvec(int *r, int **m, int *v, unsigned n1,
  unsigned n2) {
  unsigned i,j;
  for (i=0; i<n1; i++) {
    r[i] = 0;
    for (j=0; j<n2; j++)
      r[i] += (m[i][j]*v[j]);
    }
  }

/* -------------------------------------------------------------------------- */

/* multiplication of two matrices (unsigned integer version)

   m        - result matrix of dimension n1 x n3
   m1       - first matrix of dimension n1 x n2
   m2       - second matrix of dimension n2 x n3
   n1,n2,n3 - dimenstions of matrices */
void mat_umult(unsigned **m, unsigned **m1, unsigned **m2,
  unsigned n1, unsigned n2, unsigned n3) {
  unsigned i,j,k;
  for (i=0; i<n1; i++)
    for (j=0; j<n3; j++) {
      m[i][j] = 0;
      for (k=0; k<n2; k++)
        m[i][j] += (m1[i][k]*m2[k][j]);
      }
  }

/* multiplication of matrix by vector from left (unsigned version)

   r     - resulting vector of dimension n2
   m     - matrix of dimension n1 x n2
   v     - vector of dimension n1
   n1,n2 - dimensions of matrix and vectors */
void mat_umult_lvec(unsigned *r, unsigned **m, unsigned *v, unsigned n1,
  unsigned n2) {
  unsigned i,j;
  for (i=0; i<n2; i++) {
    r[i] = 0;
    for (j=0; j<n1; j++)
      r[i] += (v[j]*m[j][i]);
    }
  }

/* multiplication of matrix by vector from right (unsigned integer version)

   r     - resulting vector of dimension n1
   m     - matrix of dimension n1 x n2
   v     - vector of dimension n2
   n1,n2 - dimensions of matrix and vectors */
void mat_umult_rvec(unsigned *r, unsigned **m, unsigned *v, unsigned n1,
  unsigned n2) {
  unsigned i,j;
  for (i=0; i<n1; i++) {
    r[i] = 0;
    for (j=0; j<n2; j++)
      r[i] += (m[i][j]*v[j]);
    }
  }

/* -------------------------------------------------------------------------- */

/* multiplication of two matrices (double precision real version)

   m        - result matrix of dimension n1 x n3
   m1       - first matrix of dimension n1 x n2
   m2       - second matrix of dimension n2 x n3
   n1,n2,n3 - dimenstions of matrices */
void mat_fmult(double **m, double **m1, double **m2,
  unsigned n1, unsigned n2, unsigned n3) {
  unsigned i,j,k;
  for (i=0; i<n1; i++)
    for (j=0; j<n3; j++) {
      m[i][j] = 0.0;
      for (k=0; k<n2; k++)
        m[i][j] += (m1[i][k]*m2[k][j]);
      }
  }

/* multiplication of two matrices (complex version)

   m        - result matrix of dimension n1 x n3
   m1       - first matrix of dimension n1 x n2
   m2       - second matrix of dimension n2 x n3
   n1,n2,n3 - dimenstions of matrices */
void mat_zmult(double complex **m, double complex **m1, double complex **m2,
  unsigned n1, unsigned n2, unsigned n3) {
  unsigned i,j,k;
  for (i=0; i<n1; i++)
    for (j=0; j<n3; j++) {
      m[i][j] = 0.0;
      for (k=0; k<n2; k++)
        m[i][j] += (m1[i][k]*m2[k][j]);
      }
  }

/* multiplication of matrix by vector from left (double real version)

   r     - resulting vector of dimension n2
   m     - matrix of dimension n1 x n2
   v     - vector of dimension n1
   n1,n2 - dimensions of matrix and vectors */
void mat_fmult_lvec(double *r, double **m, double *v, unsigned n1,
  unsigned n2) {
  unsigned i,j;
  for (i=0; i<n2; i++) {
    r[i] = 0.0;
    for (j=0; j<n1; j++)
      r[i] += (v[j]*m[j][i]);
    }
  }

/* multiplication of matrix by vector from right (double real version)

   r     - resulting vector of dimension n1
   m     - matrix of dimension n1 x n2
   v     - vector of dimension n2
   n1,n2 - dimensions of matrix and vectors */
void mat_fmult_rvec(double *r, double **m, double *v, unsigned n1,
  unsigned n2) {
  unsigned i,j;
  for (i=0; i<n1; i++) {
    r[i] = 0.0;
    for (j=0; j<n2; j++)
      r[i] += (m[i][j]*v[j]);
    }
  }

/* multiply two square matrices in place, product in the first one

   a,b - two square matrices (product is saved to a)
   n   - dimension of the matrices */
void mat_fmult_sqr(double **a, double **b, unsigned n) {
  unsigned i,j,k;
  double *x;
  /* auxiliary vector */
  x = vec_falloc(n);
  /* multilication (a' = a.b) */
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      x[j] = a[i][j];
      a[i][j] = 0.0;
      for (k=0; k<n; k++) {
        if (k>j)
          a[i][j] += (a[i][k]*b[k][j]);
        else
          a[i][j] += (x[k]*b[k][j]);
        }
      }
  /* free allocated memory */
  vec_ffree(x);
  }

/* multiply two square matrices in place, product in the second one

   a,b - two square matrices (product is saved to b)
   n   - dimension of the matrices */
void mat_fmult_sql(double **a, double **b, unsigned n) {
  unsigned i,j,k;
  double *x;
  /* auxiliary vector */
  x = vec_falloc(n);
  /* multilication (a' = a.b) */
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      x[j] = b[j][i];
      b[j][i] = 0.0;
      for (k=0; k<n; k++) {
        if (k>j)
          b[j][i] += (a[j][k]*b[k][i]);
        else
          b[j][i] += (a[j][k]*x[k]);
        }
      }
  /* free allocated memory */
  vec_ffree(x);
  }

/* multiply two square 4c-format matrices in place, product in the first one

   a,b - two square matrices (product is saved to a) */
void mat_fmult_sqr_4cf(double *a, double *b) {
  static unsigned i,j,k;
  static double x[4];
  for (i=0; i<4; i++)
    for (j=0; j<4; j++) {
      x[j] = a[4*j+i];
      a[4*j+i] = 0.0;
      for (k=0; k<4; k++) {
        if (k>j)
          a[4*j+i] += (a[4*k+i]*b[4*j+k]);
        else
          a[4*j+i] += (x[k]*b[4*j+k]);
        }
      }
  }

/* multiply two square 4c-format matrices in place, product in the second one

   a,b - two square matrices (product is saved to b) */
void mat_fmult_sql_4cf(double *a, double *b) {
  static unsigned i,j,k;
  double x[4];
  for (i=0; i<4; i++)
    for (j=0; j<4; j++) {
      x[j] = b[4*i+j];
      b[4*i+j] = 0.0;
      for (k=0; k<4; k++) {
        if (k>j)
          b[4*i+j] += (a[4*k+j]*b[4*i+k]);
        else
          b[4*i+j] += (a[4*k+j]*x[k]);
        }
      }
  }

/* multiply square matrix by smaller square matrix block from the right

   a  - double precision real matrix (product will be here)
   b  - double precision real square matrix block
   na - dimension of the matrix a
   nb - dimension of the block b */
void mat_fmult_sqsub(double **a, double **b, unsigned na, unsigned nb) {
  unsigned i,j,k;
  double *x;
  /* auxiliary vector */
  x = vec_falloc(nb);
  /* multiply matrix in place */
  for (i=0; i<na; i++)
    for (j=0; j<nb; j++) {
      x[j] = a[i][j];
      a[i][j]=0.0;
      for (k=0; k<nb; k++) {
        if (k>j)
          a[i][j] += (a[i][k]*b[k][j]);
        else
          a[i][j] += (x[k]*b[k][j]);
        }
      }
  /* free memory */
  vec_ffree(x);
  }

/* similarity transformation of square matrix in place

   a - matrix which will be transformed (A'' = U^T.A.U)
   u - orthogonal matrix for transformation
   n - dimension of the matrices */
void mat_fmult_sim(double **a, double **u, unsigned n) {
  unsigned i,j,k;
  double *x;
  /* auxiliary vector */
  x = vec_falloc(n);
  /* multiply matrix from the right (A' = A.U) */
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      x[j] = a[i][j];
      a[i][j]=0.0;
      for (k=0; k<n; k++) {
        if (k>j)
          a[i][j] += (a[i][k]*u[k][j]);
        else
          a[i][j] += (x[k]*u[k][j]);
        }
      }
  /* multiply A' by U^T from the left (A'' = U^T.A') */
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      x[j] = a[j][i];
      a[j][i] = 0.0;
      for (k=0; k<n; k++) {
        if (k>j)
          a[j][i] += (u[k][j]*a[k][i]);
        else
          a[j][i] += (u[k][j]*x[k]);
        }
      }
  /* free allocated memory */
  vec_ffree(x);
  }

/* -------------------------------------------------------------------------- */
