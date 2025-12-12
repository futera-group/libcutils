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
#include "mol/distance.h"

#include <stdio.h>

/* -------------------------------------------------------------------------- */

/* wrap point coordinates to unit cell in periodic space

   x - the point coordinates (input/output)
   r - space origin
   n - space dimension
   s - shifting IDs (output)
   b - cell dimensions */
void pbc_shift(double *x, double *r, unsigned n, int *s, double *b) {
  unsigned i;
  if (b) {
    /* shift origin */
    if (r) {
      /* save shifting IDs */
      if (s) {
        for (i=0; i<n; i++) {
          if (b[i]>0.0) {
            while ((x[i]-r[i])>b[i]/2.0) {
              x[i] -= b[i];
              s[i]--;
              }
            while ((x[i]-r[i])<(-b[i]/2.0)) {
              x[i] += b[i];
              s[i]++;
              }
            }
          }
        }
      /* no shifting IDs */
      else {
        for (i=0; i<n; i++) {
          if (b[i]>0.0) {
            while ((x[i]-r[i])>b[i]/2.0) 
              x[i] -= b[i];
            while ((x[i]-r[i])<(-b[i]/2.0))
              x[i] += b[i];
            }
          }
        }
      }
    /* zero origin */
    else {
      /* save shifting IDs */
      if (s) {
        for (i=0; i<n; i++) {
          if (b[i]>0.0) {
            while (x[i]>b[i]/2.0) {
              x[i] -= b[i];
              s[i]--;
              }
            while (x[i]<(-b[i]/2.0)) {
              x[i] += b[i]; 
              s[i]++;
              }
            }
          }
        }
      /* no shifting IDs */
      else {
        for (i=0; i<n; i++) {
          if (b[i]>0.0) {
            while (x[i]>b[i]/2.0) 
              x[i] -= b[i];
            while (x[i]<(-b[i]/2.0))
              x[i] += b[i];
            }
          }
        }
      }
    }
  }

/* -------------------------------------------------------------------------- */

/* calculate distance between two points in periodic space

   c1,c2 - coordinates of the points
   n     - space dimension
   b     - periodic box dimensions */
double pbc_dist(double *c1, double *c2, unsigned n, double *b) {
  unsigned i;
  double r,r2 = 0.0;
  /* periodic boundary */
  if (b) {
    for (i=0; i<n; i++) {
      r = fabs(c1[i]-c2[i]);
      while (r>b[i]/2.0)
        r = fabs(r-b[i]);
      while (r<(-b[i]/2.0))
        r += b[i];
      r2 += r*r;
      }
    }
  /* aperiodic system */
  else
    r2 = dist_r2(c1,c2,n);
  return(sqrt(r2));
  }

/* -------------------------------------------------------------------------- */
