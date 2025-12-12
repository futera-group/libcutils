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

/* evaluate binomical function */
double math_binom_f(unsigned n, unsigned k) {
  if (k>n)
    msg_error("invalid arguments of binomical function (k>n)",1);
  if (n<171)
    return(floor(0.5+math_fact_f(n)/(math_fact_f(k)*math_fact_f(n-k))));
  return(floor(0.5+
    exp(math_fact_ln_f(n)-math_fact_ln_f(k)-math_fact_ln_f(n-k))));
  }

/* evaluate binomical function */
unsigned math_binom_u(unsigned n, unsigned k) {
  unsigned i,m,n1,n2;
  if (k>n)
    msg_error("invalid arguments of binomical function (k>n)",1);
  if (k==0)
    return(1);
  if (k==1)
    return(n);
  m = (k<(n-k) ? k : n-k);
  for (i=0,n1=1,n2=1; i<m; i++) {
    if (((n-i)%(m-i))==0)
      n1 *= (n-i)/(m-i);
    else if ((n1%(m-i))==0) {
      n1 /= (m-i);
      n1 *= (n-i);
      }
    else {
      n1 *= (n-i);
      n2 *= (m-i);
      }
    }
  return(n1/n2);
  }

/* -------------------------------------------------------------------------- */
