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
#include "cmn/message.h"

/* -------------------------------------------------------------------------- */

/* internal data for error function */
static unsigned er_ncof=28;
static double er_cof[28]={
 -1.3026537197817094E+00, 6.4196979235649026E-01,
  1.9476473204185836E-02,-9.5615147868086310E-03,
 -9.4659534448203600E-04, 3.6683949785276100E-04,
  4.2523324806907000E-05,-2.0278578112534000E-05,
 -1.6242900046470000E-06, 1.3036558355800000E-06,
  1.5626441722000000E-08,-8.5238095915000000E-08,
  6.5290544390000000E-09, 5.0593434950000000E-09,
 -9.9136415600000000E-10,-2.2736512200000000E-10,
  9.6467911000000000E-11, 2.3940380000000000E-12,
 -6.8860270000000000E-12, 8.9448700000000000E-13,
  3.1309200000000000E-13,-1.1270800000000000E-13,
  3.8100000000000000E-16, 7.1060000000000000E-15,
 -1.5230000000000000E-15,-9.4000000000000000E-17,
  1.2100000000000000E-16,-2.8000000000000000E-17};

/* -------------------------------------------------------------------------- */

/* evaluate error function */
double math_erf(double x) {
  if (x>=0.0)
    return(1.0-math_erfc_chebyshev(x));
  return(math_erfc_chebyshev(-x)-1.0);
  }

/* evaluate complementary error function */
double math_erfc(double x) {
  if (x>=0.0)
    return(math_erfc_chebyshev(x));
  else
    return(2.0-math_erfc_chebyshev(-x));
  }

/* evaluate inversion of complementary error function */
double math_erfc_inv(double p) {
  double x,err,t,pp;
  unsigned i;
  if (p>=2.0)
    return(-100.0);
  if (p<=0.0)
    return(100.0);
  pp = (p<1.0 ? p : 2.0-p);
  t = sqrt(-2.0*log(pp/2.0));
  x =- 0.70711*((2.30753+t*0.27061)/(1.0+t*(0.99229+t*0.04481))-t);
  for (i=0; i<2; i++) {
    err = math_erfc(x)-pp;
    x += err/(1.12837916709551257*exp(-(x*x))-x*err);
    }
  return(p<1.0 ? x : -x);
  }

/* evaluate complementary error function by Chebyshev interpolation */
double math_erfc_chebyshev(double z) {
  double t,ty,tmp,d = 0.0,dd = 0.0;
  int i;
  if (z<0.0)
    msg_error("negative argument in complementary error function",1);
  t = 2.0/(2.0+z);
  ty = 4.0*t-2.0;
  for (i=er_ncof-1; i>0; i--) {
    tmp = d;
    d = ty*d-dd+er_cof[i];
    dd = tmp;
    }
  return(t*exp(-z*z+0.5*(er_cof[0]+ty*d)-dd));
  }

/* -------------------------------------------------------------------------- */
