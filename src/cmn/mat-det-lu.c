/******************************************************************************\
 *                                                                            * 
 *  Libcutils - library of C function                                         * 
 *                                                                            *
 *  Version:             3.4                                                  * 
 *  Date:                28/01/2017                                           *
 *                                                                            * 
 *  Author:              Zdenek Futera                                        * 
 *                                                                            * 
 *  Address: University of South Bohemia                                      *
 *           Faculty of Science, Institute of Physics                         *
 *           Branisovska 1760, 370 05 Ceske Budejovice                        *
 *           Czech Republic                                                   *
 *                                                                            *
 *  Email:   zfutera@prf.jcu.cz                                               *
 *                                                                            * 
\******************************************************************************/

#include <complex.h>

/* -------------------------------------------------------------------------- */

/* calculate matrix determinant by LU decomposition (double-precision real)

   m - the LU-decomposed matrix
   s - number of exchanged rows in LU decomposition
   n - dimension of the matrix */
double mat_fdet_lu(double **m, unsigned s, unsigned n) {
  double d = 1.0;
  unsigned i;
  /* diagonal elements */
  for (i=0; i<n; i++)
    d *= m[i][i];
  /* determinant */
  return(s%2 ? -d : d);
  }

/* calculate matrix determinant by LU decomposition (double-precision complex)

   m - the LU-decomposed matrix
   s - number of exchanged rows in LU decomposition
   n - dimension of the matrix */
complex double mat_zdet_lu(complex double **m, unsigned s, unsigned n) {
  complex double d = 1.0;
  unsigned i;
  /* diagonal elements */
  for (i=0; i<n; i++)
    d *= m[i][i];
  /* determinant */
  return(s%2 ? -d : d);
  }

/* -------------------------------------------------------------------------- */
