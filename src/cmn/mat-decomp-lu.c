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
#include <math.h>
#include <cmn/message.h>

/* -------------------------------------------------------------------------- */

/* LU decomposition of a square matrix (double-precision real)
 
   m - the matrix
   n - size of the matrix
   p - permutation index array (output)
   s - number of exchanged rows (output) */
void mat_fdecomp_lu(double **m, unsigned n, unsigned *p, unsigned *s) {
  unsigned i,j,k,ix;
  double mx,t,*q;
  /* initialization */
  for (i=0; i<n; i++)
    p[i] = i;
  (*s) = 0;
  /* decomposition */
  for (i=0; i<n; i++) {
    /* scaling factor */
    mx = 0.0;
    ix = i;
    for (j=i; j<n; j++) {
      t = fabs(m[j][i]);
      if (t > mx) {
        mx = t;
        ix = j;
        }
      }
    /* sanity check */
    if (mx == 0.0)
      msg_error("singular matrix in LU decomposition",1);
    /* pivoting */
    if (ix != i) {
      /* row IDs */
      j = p[i];
      p[i] = p[ix];
      p[ix] = j;
      /* row values */
      q = m[i];
      m[i] = m[ix];
      m[ix] = q;
      /* counter */
      (*s)++;
      }
    /* elimination */
    for (j=i+1; j<n; j++) {
      m[j][i] /= m[i][i];
      for (k=i+1; k<n; k++)
        m[j][k] -= m[j][i] * m[i][k];
      }
    }
  }

/* LU decomposition of a square matrix (double-precision complex)
 
   m - the matrix
   n - size of the matrix
   p - permutation index array (output)
   s - number of exchanged rows (output) */
void mat_zdecomp_lu(complex double **m, unsigned n, unsigned *p, unsigned *s) {
  unsigned i,j,k,ix;
  double mx,t;
  complex double *q;
  /* initialization */
  for (i=0; i<n; i++)
    p[i] = i;
  (*s) = 0;
  /* decomposition */
  for (i=0; i<n; i++) {
    /* scaling factor */
    mx = 0.0;
    ix = i;
    for (j=i; j<n; j++) {
      t = cabs(m[j][i]);
      if (t > mx) {
        mx = t;
        ix = j;
        }
      }
    /* sanity check */
    if (mx == 0.0)
      msg_error("singular matrix in LU decomposition",1);
    /* pivoting */
    if (ix != i) {
      /* row IDs */
      j = p[i];
      p[i] = p[ix];
      p[ix] = j;
      /* row values */
      q = m[i];
      m[i] = m[ix];
      m[ix] = q;
      /* counter */
      (*s)++;
      }
    /* elimination */
    for (j=i+1; j<n; j++) {
      m[j][i] /= m[i][i];
      for (k=i+1; k<n; k++)
        m[j][k] -= m[j][i] * m[i][k];
      }
    }
  }

/* -------------------------------------------------------------------------- */
