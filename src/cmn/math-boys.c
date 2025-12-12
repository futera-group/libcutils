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

#define BOYS_EPS   1.0E-10  /* minimum threshold in Boys function */

/* prefactor of Boys function f(n,x): (2n-1)!!*sqrt(pi)/(2^(n+1)) */
#define BOYS_PREF_00      0.88622692545275805198 /* prefactor for 0th order */
#define BOYS_PREF_01      0.44311346272637902599 /* prefactor for 1st order */
#define BOYS_PREF_02      0.66467019408956851123 /* prefactor for 2nd order */
#define BOYS_PREF_03      1.66167548522392127808 /* prefactor for 3rd order */
#define BOYS_PREF_04      5.81586419828372402918 /* prefactor for 4th order */
#define BOYS_PREF_05     26.17138889227675946358 /* prefactor for 5th order */
#define BOYS_PREF_06    143.94263890752216639157 /* prefactor for 6th order */
#define BOYS_PREF_07    935.62715289889422365377 /* prefactor for 7th order */
#define BOYS_PREF_08   7017.20364674170650687302 /* prefactor for 8th order */
#define BOYS_PREF_09  59646.23099730450485367328 /* prefactor for 9th order */
#define BOYS_PREF_10 566639.19447439280338585377 /* prefactor for 10th order */

/* -------------------------------------------------------------------------- */

/* evaluate integral Boys function (core of GTO integrals) */
double math_boys(unsigned n, double t) {
  if (t<0.0)
    msg_error("negative argument in Boys function",1);
  if (n>10) {
    if (t<BOYS_EPS)
      return(1.0/(2.0*n+1));
    return(((2.0*n-1.0)*math_boys(n-1,t)-exp(-t))/(2.0*t));
    }
  switch (n) {
    case 0:  /* zero-th order */
      if (t<BOYS_EPS)
        return(1.0);
      return(BOYS_PREF_00/pow(t,0.5)*erf(sqrt(t)));
    case 1:  /* first order */
      if (t<BOYS_EPS)
        return(MATH_INV_03);
      return(BOYS_PREF_01/pow(t,1.5)*math_igamma_p(1.5,t));
    case 2:  /* second order */
      if (t<BOYS_EPS)
        return(0.2);
      return(BOYS_PREF_02/pow(t,2.5)*math_igamma_p(2.5,t));
    case 3:  /* third order */
      if (t<BOYS_EPS)
        return(MATH_INV_07);
      return(BOYS_PREF_03/pow(t,3.5)*math_igamma_p(3.5,t));
    case 4:  /* fourth order */
      if (t<BOYS_EPS)
        return(MATH_INV_09);
      return(BOYS_PREF_04/pow(t,4.5)*math_igamma_p(4.5,t));
    case 5:  /* fifth order */
      if (t<BOYS_EPS)
        return(MATH_INV_11);
      return(BOYS_PREF_05/pow(t,5.5)*math_igamma_p(5.5,t));
    case 6:  /* sixth order */
      if (t<BOYS_EPS)
        return(MATH_INV_13);
      return(BOYS_PREF_06/pow(t,6.5)*math_igamma_p(6.5,t));
    case 7:  /* seventh order */
      if (t<BOYS_EPS)
        return(MATH_INV_15);
      return(BOYS_PREF_07/pow(t,7.5)*math_igamma_p(7.5,t));
    case 8:  /* eith order */
      if (t<BOYS_EPS)
        return(MATH_INV_17);
      return(BOYS_PREF_08/pow(t,8.5)*math_igamma_p(8.5,t));
    case 9:  /* nine order */
      if (t<BOYS_EPS)
        return(MATH_INV_19);
      return(BOYS_PREF_09/pow(t,9.5)*math_igamma_p(9.5,t));
    case 10: /* tenth order */
      if (t<BOYS_EPS)
        return(MATH_INV_21);
      return(BOYS_PREF_10/pow(t,10.5)*math_igamma_p(10.5,t));
    }  
  return(0.0);
  }

/* -------------------------------------------------------------------------- */
