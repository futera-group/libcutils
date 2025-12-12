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

/* -------------------------------------------------------------------------- */

/* calculate arithmetic average and variance from array of values
 
   x - the array with values
   n - number of values in the array
   s - calculated variance */
double math_avrg(double *x, long unsigned n, double *s) {
  double avrg=0.0,sgm=0.0;
  long unsigned i;
  for (i=0; i<n; i++)
    avrg += x[i];
  avrg /= n;
  if (s) {
    for (i=0; i<n; i++)
      sgm += pow(x[i]-avrg,2.0);
    if (n>1)
      sgm = sqrt(sgm/(n-1));
    (*s) = sgm;
    }
  return(avrg);
  }

/* -------------------------------------------------------------------------- */
