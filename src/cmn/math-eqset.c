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
#include "cmn/message.h"

/* -------------------------------------------------------------------------- */

#define EQSET_ZERO 1.0E-10

/* -------------------------------------------------------------------------- */

/* change rows in matrix so that the biggest element in the first column
   is in the first row

   mat - matrix representing the equation set
   n   - number of equations
   id  - submatrix ID */
void math_eqset_gauss_rows(double **mat, unsigned n, unsigned id) {
  unsigned i,mrow=id;
  double val=0.0;
  for (i=id; i<n; i++)
    if (fabs(mat[id][i])>val) { 
      mrow = i;
      val = fabs(mat[id][i]);
      }
  if (mrow>id) {
    for (i=id; i<n+1; i++) {
      val = mat[id][i];
      mat[id][i] = mat[mrow][i];
      mat[mrow][i] = val;
      }
    }
  }

/* recursive Gauss elimination

   mat - matrix representing the equation set
   n   - number of equations
   id  - submatrix ID */
void math_eqset_gauss_elim(double **mat, unsigned n, unsigned id) {
  unsigned i,j;
  double v;
  if (id==n-1) {
    if (fabs(mat[id][id])<EQSET_ZERO) {
      if (fabs(mat[id][id+1])<EQSET_ZERO)
        msg_error("reduction of equations in gauss elimination",1);
      else {
        msg_error("set of linear equations has no solution",1);
        }
      }
    mat[id][n] /= mat[id][id];
    mat[id][id] = 1.0;
    }
  else {
    math_eqset_gauss_rows(mat,n,id);   
    if (fabs(mat[id][id])<EQSET_ZERO) 
      msg_error("reduction of variables in gauss elimination",1);
    for (i=id+1; i<n+1; i++)
      mat[id][i] /= mat[id][id];
    mat[id][id] = 1.0;
    for (i=id+1; i<n; i++) {
      v = mat[i][id];
      for (j=id; j<n+1; j++) 
        mat[i][j] = mat[i][j]-mat[id][j]*v;
      }
    math_eqset_gauss_elim(mat,n,id+1);
    }
  }

/* solve set of linear equations by Gauss elimination

   mat - matrix representing the equation set
   x   - vector for the solution
   n   - number of equations */
void math_eqset_gauss(double **mat, double *x, unsigned n) {
  int i;
  unsigned j;
  math_eqset_gauss_elim(mat,n,0);
  for (i=n-1; i>=0; i--) {
    x[i] = mat[i][n];
    for (j=i+1; j<n; j++)
      x[i] -= mat[i][j]*x[j];
      }
  }

/* -------------------------------------------------------------------------- */
