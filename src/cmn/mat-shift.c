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

/* shift rows down in matrix (integer)

   m     - pointer to integer matrix
   n1,n2 - dimension of the matrix */
void mat_idshift(int **m, unsigned n1, unsigned n2) {
  long int i,j;
  for (i=n1-1; i>0; i--)
    for (j=0; j<n2; j++)
      m[i][j] = m[i-1][j];
  for (j=0; j<n2; j++)
    m[0][j] = 0;
  }

/* shift rows down in matrix (unsigned integer)

   m     - pointer to unsigned integer matrix
   n1,n2 - dimension of the matrix */
void mat_udshift(unsigned **m, unsigned n1, unsigned n2) {
  long int i,j;
  for (i=n1-1; i>0; i--)
    for (j=0; j<n2; j++)
      m[i][j] = m[i-1][j];
  for (j=0; j<n2; j++)
    m[0][j] = 0;
  }

/* shift rows down in matrix (double precision real)

   m     - pointer to double precision real matrix
   n1,n2 - dimension of the matrix */
void mat_fdshift(double **m, unsigned n1, unsigned n2) {
  long int i,j;
  for (i=n1-1; i>0; i--)
    for (j=0; j<n2; j++)
      m[i][j] = m[i-1][j];
  for (j=0; j<n2; j++)
    m[0][j] = 0.0;
  }

/* -------------------------------------------------------------------------- */

/* shift rows up in matrix (integer)

   m     - pointer to integer matrix
   n1,n2 - dimension of the matrix */
void mat_iushift(int **m, unsigned n1, unsigned n2) {
  long int i,j;
  for (i=0; i<n1-1; i++)
    for (j=0; j<n2; j++)
      m[i][j] = m[i+1][j];
  for (j=0; j<n2; j++)
    m[n1-1][j] = 0;
  }

/* shift rows up in matrix (unsigned integer)

   m     - pointer to unsigned integer matrix
   n1,n2 - dimension of the matrix */
void mat_uushift(unsigned **m, unsigned n1, unsigned n2) {
  long int i,j;
  for (i=0; i<n1-1; i++)
    for (j=0; j<n2; j++)
      m[i][j] = m[i+1][j];
  for (j=0; j<n2; j++)
    m[n1-1][j] = 0;
  }

/* shift rows up in matrix (double precision real)

   m     - pointer to double precision real matrix
   n1,n2 - dimension of the matrix */
void mat_fushift(double **m, unsigned n1, unsigned n2) {
  long int i,j;
  for (i=0; i<n1-1; i++)
    for (j=0; j<n2; j++)
      m[i][j] = m[i+1][j];
  for (j=0; j<n2; j++)
    m[n1-1][j] = 0.0;
  }

/* -------------------------------------------------------------------------- */

/* shift columns left in matrix (integer)

   m     - pointer to integer real matrix
   n1,n2 - dimension of the matrix */
void mat_ilshift(int **m, unsigned n1, unsigned n2) {
  long int i,j;
  for (i=0; i<n1; i++)
    for (j=0; j<n2-1; j--)
      m[i][j] = m[i][j+1];
  for (i=0; i<n1; i++)
    m[i][n2-1] = 0;
  }

/* shift columns left in matrix (unsigned integer)

   m     - pointer to unsigned integer real matrix
   n1,n2 - dimension of the matrix */
void mat_ulshift(unsigned **m, unsigned n1, unsigned n2) {
  long int i,j;
  for (i=0; i<n1; i++)
    for (j=0; j<n2-1; j--)
      m[i][j] = m[i][j+1];
  for (i=0; i<n1; i++)
    m[i][n2-1] = 0;
  }

/* shift columns left in matrix (double precision real)

   m     - pointer to double precision real matrix
   n1,n2 - dimension of the matrix */
void mat_flshift(double **m, unsigned n1, unsigned n2) {
  long int i,j;
  for (i=0; i<n1; i++)
    for (j=0; j<n2-1; j--)
      m[i][j] = m[i][j+1];
  for (i=0; i<n1; i++)
    m[i][n2-1] = 0.0;
  }

/* -------------------------------------------------------------------------- */

/* shift columns right in matrix (integer)

   m     - pointer to integer matrix
   n1,n2 - dimension of the matrix */
void mat_irshift(int **m, unsigned n1, unsigned n2) {
  long int i,j;
  for (i=0; i<n1; i++)
    for (j=n2-1; j>0; j--)
      m[i][j] = m[i][j-1];
  for (i=0; i<n1; i++)
    m[i][0] = 0;
  }

/* shift columns right in matrix (unsigned integer)

   m     - pointer to unsigned integer matrix
   n1,n2 - dimension of the matrix */
void mat_urshift(unsigned **m, unsigned n1, unsigned n2) {
  long int i,j;
  for (i=0; i<n1; i++)
    for (j=n2-1; j>0; j--)
      m[i][j] = m[i][j-1];
  for (i=0; i<n1; i++)
    m[i][0] = 0;
  }

/* shift columns right in matrix (double precision real)

   m     - pointer to double precision real matrix
   n1,n2 - dimension of the matrix */
void mat_frshift(double **m, unsigned n1, unsigned n2) {
  long int i,j;
  for (i=0; i<n1; i++)
    for (j=n2-1; j>0; j--)
      m[i][j] = m[i][j-1];
  for (i=0; i<n1; i++)
    m[i][0] = 0.0;
  }

/* -------------------------------------------------------------------------- */
