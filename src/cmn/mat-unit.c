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

/* set square matrix to unit (int version)

   m - the real square matrix
   n - diminsion of the matrix */
void mat_iunit(int **m, unsigned n) {
  unsigned i,j;
  for (i=0; i<n; i++) {
    m[i][i] = 1;
    for (j=i+1; j<n; j++)
      m[i][j] = m[j][i] = 0;
    }
  }

/* set square matrix to unit (unsigned int version)

   m - the real square matrix
   n - diminsion of the matrix */
void mat_uunit(unsigned **m, unsigned n) {
  unsigned i,j;
  for (i=0; i<n; i++) {
    m[i][i] = 1;
    for (j=i+1; j<n; j++)
      m[i][j] = m[j][i] = 0;
    }
  }

/* set square matrix to unit (double precision real version)

   m - the real square matrix
   n - diminsion of the matrix */
void mat_funit(double **m, unsigned n) {
  unsigned i,j;
  for (i=0; i<n; i++) {
    m[i][i] = 1.0;
    for (j=i+1; j<n; j++)
      m[i][j] = m[j][i] = 0.0;
    }
  }

/* -------------------------------------------------------------------------- */
