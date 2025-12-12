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

#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* calculate histogram from given values

   x  - vector with values
   n  - number of values
   mn - first bin
   mx - last bin
   nb - number of bins
   s  - step (bin width) */
unsigned long *math_hist_calc(double *x, long unsigned n, double mn, double mx,
  unsigned nb, double *s) {
  long unsigned i,*h;
  double st;
  h = vec_lualloc(nb);
  vec_luset(h,0,nb);
  st = (mx-mn)/(nb-1);
  for (i=0; i<n; i++)
    h[(unsigned)((x[i]-mn)/st)]++;
  if (s)
    (*s) = st;
  return(h);
  }

/* create normalized histogram from number of counts

   h - array with unnormalized histogram
   n - number of bins in the histogram
   d - bin width */
double *math_hist_norm(unsigned long *h, unsigned n, double d) {
  double *g=NULL,s=0.0;
  unsigned i;
  g = vec_falloc(n);
  for (i=0; i<n-1; i++)
    s += (0.5*d*(h[i]+h[i+1]));
  for (i=0; i<n; i++)
    g[i] = ((double)h[i])/s;
  return(g);
  }

/* -------------------------------------------------------------------------- */
