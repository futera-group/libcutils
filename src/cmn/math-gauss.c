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
#include "cmn/math.h"

#define MATH_GAUSS_MAX_ITER 10000
#define MATH_GAUSS_CNV_CRIT 1.0E-7

/* -------------------------------------------------------------------------- */

/* evaluate normalized gaussian function
 
   m - mean value
   s - sigma
   x - functional value */
double math_gauss_fce(double m, double s, double x) {
  return(exp(-0.5*(x-m)*(x-m)/(s*s))/(s*sqrt(2.0*M_PI)));
  }

/* fit values to normalized gaussian function

   x - array with X values
   y - array with Y values
   n - number of values
   m - fitted mean value
   s - fitted half-width */
void math_gauss_fit(double *x, double *y, unsigned n, double *m, double *s) {
  double rr,r0,c11,c12,c22,b1,b2,ds,dm,db,gm,gs,gx;
  double mu,sigma,min_m,min_s;
  unsigned i,it;
  /* initialization */
  rr = 0.0;
  mu = math_avrg(x,n,&sigma);
  min_s = sigma;
  min_m = mu;
  /* gaussian fit */
  for (it=0; it<MATH_GAUSS_MAX_ITER; it++) {
    r0 = rr;
    c11 = c12 = c22 = b1 = b2 = rr = 0.0;
    for (i=0; i<n; i++) {
      gx = math_gauss_fce(mu,sigma,x[i]);
      gm = (x[i]-mu)*gx/pow(sigma,2.0);
      gs = pow(x[i]-mu,2.0)*gx/pow(sigma,3.0);
      db = y[i]-gx;
      rr += (db*db);
      b1 += (gm*db);
      b2 += (gs*db);
      c11 += (gm*gm);
      c12 += (gm*gs);
      c22 += (gs*gs);
      }
    /* mu and sigma correction */
    ds = (b2*c11-b1*c12)/(c11*c22-c12*c12);
    dm = (b1-ds*c12)/c11;
    mu += dm;
    sigma += ds;
    /* check convergency */
    if (it) {
      if (rr<r0) {
        min_m = mu;
        min_s = sigma;
        }
      if (fabs(r0-rr)<MATH_GAUSS_CNV_CRIT)
        break;
      }
    }
  /* return best fit */
  (*m) = min_m;
  (*s) = min_s;
  }

/* -------------------------------------------------------------------------- */
