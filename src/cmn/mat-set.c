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

/* set all elements of matrix to one value (integer version)

   mat  - pointer to the matrix
   nrow - number of rows in the matrix
   ncol - number of columns in the matrix
   val  - the value for setting */
void mat_iset(int **mat, int val, unsigned nrow, unsigned ncol) {
  unsigned i,j;
  for (i=0; i<nrow; i++)
    for (j=0; j<ncol; j++)
      mat[i][j] = val;
  }

/* set all elements of matrix to one value (unsigned int version)

   mat  - pointer to the matrix
   nrow - number of rows in the matrix
   ncol - number of columns in the matrix
   val  - the value for setting */
void mat_uset(unsigned **mat, unsigned val, unsigned nrow, unsigned ncol) {
  unsigned i,j;
  for (i=0; i<nrow; i++)
    for (j=0; j<ncol; j++)
      mat[i][j] = val;
  }

/* set all elements of matrix to one value (short integer version)

   mat  - pointer to the matrix
   nrow - number of rows in the matrix
   ncol - number of columns in the matrix
   val  - the value for setting */
void mat_siset(short int **mat, short int val, unsigned nrow, unsigned ncol) {
  unsigned i,j;
  for (i=0; i<nrow; i++)
    for (j=0; j<ncol; j++)
      mat[i][j] = val;
  }

/* set all elements of matrix to one value (unsigned long int version)

   mat  - pointer to the matrix
   nrow - number of rows in the matrix
   ncol - number of columns in the matrix
   val  - the value for setting */
void mat_luset(long unsigned **mat, long unsigned val,
  unsigned nrow, unsigned ncol) {
  unsigned i,j;
  for (i=0; i<nrow; i++)
    for (j=0; j<ncol; j++)
      mat[i][j] = val;
  }

/* set all elements of matrix to one value (double precision real version)

   mat  - pointer to the matrix
   nrow - number of rows in the matrix
   ncol - number of columns in the matrix
   val  - the value for setting */
void mat_fset(double **mat, double val, unsigned nrow, unsigned ncol) {
  unsigned i,j;
  for (i=0; i<nrow; i++)
    for (j=0; j<ncol; j++)
      mat[i][j] = val;
  }

/* set all elements of matrix to one value (complex version)

   m  - pointer to the matrix
   nr - number of rows in the matrix
   nc - number of columns in the matrix
   v  - the value for setting */
void mat_zset(double complex **m, double complex v, unsigned nr, unsigned nc) {
  unsigned i,j;
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
      m[i][j] = v;
  }

/* -------------------------------------------------------------------------- */
