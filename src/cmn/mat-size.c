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

#include <stdlib.h>
#include "cmn/matrix.h"

/* -------------------------------------------------------------------------- */

/* change number of rows in unsigned integer matrix

   m  - pointer to the matrix
   c  - number of columns
   r1 - number of rows before resizing
   r2 - number of rows after resizing */
unsigned **mat_uresize_r(unsigned **m, unsigned c, unsigned r1, unsigned r2) {
  unsigned i,j,n,**x = NULL;
  if (r1==r2)
    x = m;
  else if (r2!=0) {
    x = mat_ualloc(r2,c);
    n = (r1<r2 ? r1 : r2);
    for (i=0; i<n; i++)
      for (j=0; j<c; j++)
        x[i][j] = m[i][j];
    mat_ufree(m,r1);
    }
  else 
    mat_ufree(m,r1);
  return(x);
  }

/* change number of rows in floating-point matrix

   m  - pointer to the matrix
   c  - number of columns
   r1 - number of rows before resizing
   r2 - number of rows after resizing */
double **mat_fresize_r(double **m, unsigned c, unsigned r1, unsigned r2) {
  double **x = NULL;
  unsigned i,j,n;
  if (r1==r2)
    x = m;
  else if (r2!=0) {
    x = mat_falloc(r2,c);
    n = (r1<r2 ? r1 : r2);
    for (i=0; i<n; i++)
      for (j=0; j<c; j++)
        x[i][j] = m[i][j];
    mat_ffree(m,r1);
    }
  else 
    mat_ffree(m,r1);
  return(x);
  }

/* -------------------------------------------------------------------------- */

/* delete specific row in unsigned integer matrix

   m - pointer to the matrix
   c - number of columns
   r - total number of rows
   d - the row for deletion */
unsigned **mat_udelete_r(unsigned **m, unsigned c, unsigned r, unsigned d) {
  unsigned i,j,**x = NULL;
  if (r>1) {
    x = mat_ualloc(r-1,c);
    for (i=0; i<d; i++)
      for (j=0; j<c; j++)
        x[i][j] = m[i][j];
    for (i=d; i<r-1; i++)
      for (j=0; j<c; j++)
        x[i][j] = m[i+1][j];
    mat_ufree(m,r);
    }
  else
    mat_ufree(m,r);
  return(x);
  }

/* delete specific row in floating-point matrix

   m - pointer to the matrix
   c - number of columns
   r - total number of rows
   d - the row for deletion */
double **mat_fdelete_r(double **m, unsigned c, unsigned r, unsigned d) {
  double **x = NULL;
  unsigned i,j;
  if (r>1) {
    x = mat_falloc(r-1,c);
    for (i=0; i<d; i++)
      for (j=0; j<c; j++)
        x[i][j] = m[i][j];
    for (i=d; i<r-1; i++)
      for (j=0; j<c; j++)
        x[i][j] = m[i+1][j];
    mat_ffree(m,r);
    }
  else
    mat_ffree(m,r);
  return(x);
  }

/* -------------------------------------------------------------------------- */
