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
#include <stdlib.h>
#include "cmn/matrix.h"
#include "cmn/message.h"

/* -------------------------------------------------------------------------- */

/* copy values from one matrix to another (integer version)

   m1 - matrix where the values are writen
   m2 - matrix where the values are read from 
   n   - number of rows in the matrices
   m   - number of columns in the matrices */
void mat_icopy(int **m1, int **m2, unsigned n, unsigned m) {
  unsigned i,j;
  if (!m1 || !m2)
    msg_error("null pointer in copy int matrix",1);
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      m1[i][j] = m2[i][j];
  }

/* create new matrix and fill it with values from another one (int version)

   mat - matrix where the values are read from 
   n   - number of rows in the matrix
   m   - number of columns in the matrix */
int **mat_icopy_new(int **mat, unsigned n, unsigned m) {
  int **nm;
  unsigned i,j;
  if (!mat)
    return(NULL);
  nm = mat_ialloc(n,m);
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      nm[i][j] = mat[i][j];
  return(nm);
  }

/* -------------------------------------------------------------------------- */

/* copy values from one matrix to another (unsigned integer version)

   m1 - matrix where the values are writen
   m2 - matrix where the values are read from 
   n   - number of rows in the matrices
   m   - number of columns in the matrices */
void mat_ucopy(unsigned **m1, unsigned **m2, unsigned n, unsigned m) {
  unsigned i,j;
  if (!m1 || !m2)
    msg_error("null pointer in copy unsigned int matrix",1);
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      m1[i][j] = m2[i][j];
  }

/* create new matrix and fill it with values from another one (uint version)

   mat - matrix where the values are read from 
   n   - number of rows in the matrix
   m   - number of columns in the matrix */
unsigned **mat_ucopy_new(unsigned **mat, unsigned n, unsigned m) {
  unsigned **nm;
  unsigned i,j;
  if (!mat)
    return(NULL);
  nm = mat_ualloc(n,m);
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      nm[i][j] = mat[i][j];
  return(nm);
  }

/* -------------------------------------------------------------------------- */

/* copy values from one matrix to another (short integer version)

   m1 - matrix where the values are writen
   m2 - matrix where the values are read from 
   n   - number of rows in the matrices
   m   - number of columns in the matrices */
void mat_sicopy(short int **m1, short int **m2, unsigned n, unsigned m) {
  unsigned i,j;
  if (!m1 || !m2)
    msg_error("null pointer in copy int matrix",1);
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      m1[i][j] = m2[i][j];
  }

/* create new matrix and fill it with values from another one (s-int version)

   mat - matrix where the values are read from 
   n   - number of rows in the matrix
   m   - number of columns in the matrix */
short int **mat_sicopy_new(short int **mat, unsigned n, unsigned m) {
  short int **nm;
  unsigned i,j;
  if (!mat)
    return(NULL);
  nm = mat_sialloc(n,m);
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      nm[i][j] = mat[i][j];
  return(nm);
  }

/* -------------------------------------------------------------------------- */

/* copy values from one matrix to another (double real version)

   m1 - matrix where the values are writen
   m2 - matrix where the values are read from 
   n   - number of rows in the matrices
   m   - number of columns in the matrices */
void mat_fcopy(double **m1, double **m2, unsigned n, unsigned m) {
  unsigned i,j;
  if (!m1 || !m2)
    msg_error("null pointer in copy double real matrix",1);
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      m1[i][j] = m2[i][j];
  }

/* copy transposed values from one matrix to another (double real version)

   m1 - matrix where the values are writen
   m2 - matrix where the values are read from 
   n   - number of rows in the matrices
   m   - number of columns in the matrices */
void mat_fcopy_trans(double **m1, double **m2, unsigned n, unsigned m) {
  unsigned i,j;
  if (!m1 || !m2)
    msg_error("null pointer in copy double real matrix",1);
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      m1[i][j] = m2[j][i];
  }

/* create new matrix and fill it with values from another one (double version)

   mat - matrix where the values are read from 
   n   - number of rows in the matrix
   m   - number of columns in the matrix */
double **mat_fcopy_new(double **mat, unsigned n, unsigned m) {
  double **nm;
  unsigned i,j;
  if (!mat)
    return(NULL);
  nm = mat_falloc(n,m);
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      nm[i][j] = mat[i][j];
  return(nm);
  }

/* create new matrix and fill it with transposed values from another one (double version)

   mat - matrix where the values are read from 
   n   - number of rows in the matrix
   m   - number of columns in the matrix */
double **mat_fcopy_trans_new(double **mat, unsigned n, unsigned m) {
  double **nm;
  unsigned i,j;
  if (!mat)
    return(NULL);
  nm = mat_falloc(m,n);
  for (i=0; i<m; i++)
    for (j=0; j<n; j++)
      nm[i][j] = mat[j][i];
  return(nm);
  }

/* copy values from sub-matrix (double real version)

   m1 - matrix where the values are writen (dimension of submatrix)
   m2 - matrix where the values are read from
   r1,c1 - first row and column in the submatrix
   r2,c2 - last row and column in the submatrix */
void mat_fcopy_sub(double **m1, double **m2, unsigned r1, unsigned c1,
  unsigned r2, unsigned c2) {
  unsigned i,j,k,l;
  if (!m1 || !m2)
    msg_error("null pointer in copy double real matrix",1);
  k = l = 0;
  for (i=r1; i<r2; i++, k++) {
    for (j=c1; j<c2; j++, l++) 
      m1[k][l] = m2[i][j];
    l = 0;
    }
  }

/* -------------------------------------------------------------------------- */

/* copy values from one matrix to another (complex real version)

   m1 - matrix where the values are writen
   m2 - matrix where the values are read from 
   n   - number of rows in the matrices
   m   - number of columns in the matrices */
void mat_zcopy(complex double **m1, complex double **m2,
  unsigned n, unsigned m) {
  unsigned i,j;
  if (!m1 || !m2)
    msg_error("null pointer in copy complex real matrix",1);
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      m1[i][j] = m2[i][j];
  }

/* -------------------------------------------------------------------------- */
