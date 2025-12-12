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

/* calculate length of hypotenuse by pythagorean theorem */
double math_fpyth(double a, double b) {
  double fa,fb;
  fa = fabs(a);
  fb = fabs(b);
  if (fa>fb)
    return(fa*sqrt(1.0+pow(fb/fa,2.0)));
  if (!fb)
    return(0.0);
  return(fb*sqrt(1.0+pow(fa/fb,2.0)));
  }

/* -------------------------------------------------------------------------- */
