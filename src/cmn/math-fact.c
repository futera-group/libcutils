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

#include "cmn/math.h"
#include "cmn/message.h"

/* -------------------------------------------------------------------------- */

/* evaluate factorial function - singed integer version */
double math_fact_if(int n) {
  if (n>0)
    return(math_fact_f(n));
  return(1.0);
  }

/* evaluate factorial function */
double math_fact_f(unsigned n) {
  static double v[171];
  static short init=1;
  unsigned i;
  if (init) {
    init = 0;
    v[0] = 1.0;
    for (i=1; i<171; i++)
      v[i] = i*v[i-1];
    }
  if (n>170)
    msg_error("argument of factorial function out of range",1);
  return v[n];
  }

/* evaluate odd factorial function */
double math_fact_odd_f(unsigned n) {
  static double v[171];
  static short init=1;
  unsigned i;
  if (!(n%2))
    msg_error("even number passed to odd factorial function",1);
  if (init) {
    init=0;
    v[0] = v[1] = v[2] = 1.0;
    for (i=3; i<171; i++) {
      if (i%2)
        v[i] = i*v[i-2];
      else
        v[i] = 1.0;
      }
    }
  if (n>170)
    msg_error("argument of factorial function out of range",1);
  return v[n];
  }

/* evaluate logarithm of factorial function */
double math_fact_ln_f(unsigned n) {
  static double v[5000];
  static short init = 1;
  unsigned i;
  if (init) {
    init = 0;
    for (i=0; i<5000; i++)
      v[i] = math_gamma_ln(i+1.0);
    }
  if (n<5000)
    return(v[n]);
  return(math_gamma_ln(n+1.0));
  }

/* -------------------------------------------------------------------------- */
