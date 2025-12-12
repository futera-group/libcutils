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

/* -------------------------------------------------------------------------- */

/* transpose integer square matrix

   mat - integer square matrix
   n   - dimension of the matrix */
void mat_itrans(int **mat, unsigned n) {
  int t;
  unsigned i,j;
  if (n==1)
    return;
  for (i=0; i<n-1; i++)
    for (j=i+1; j<n; j++) {
      t = mat[i][j];
      mat[i][j] = mat[j][i];
      mat[j][i] = t;
      }
  }

/* transpose unsigned integer square matrix

   mat - unsigned integer square matrix
   n   - dimension of the matrix */
void mat_utrans(unsigned **mat, unsigned n) {
  unsigned i,j,t;
  if (n==1)
    return;
  for (i=0; i<n-1; i++)
    for (j=i+1; j<n; j++) {
      t = mat[i][j];
      mat[i][j] = mat[j][i];
      mat[j][i] = t;
      }
  }

/* transpose real square matrix

   mat - double precision real square matrix
   n   - dimension of the matrix */
void mat_ftrans(double **mat, unsigned n) {
  unsigned i,j;
  double t;
  if (n==1)
    return;
  for (i=0; i<n-1; i++)
    for (j=i+1; j<n; j++) {
      t = mat[i][j];
      mat[i][j] = mat[j][i];
      mat[j][i] = t;
      }
  }

/* -------------------------------------------------------------------------- */
