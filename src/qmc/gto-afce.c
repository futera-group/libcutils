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
#include <cmn/math.h>
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* evaluate auxiliary f(l,m,a,b) function

   id  - integer specifier
   n   - order of the function
   l,m - specific Cartesian components of GTOs
   a,b - specific Cartesian coordinate of GTOs */
double gto_afce_f(unsigned n, unsigned l, unsigned m, double a, double b) {
  unsigned i,mn,mx;
  double b1,b2,ss;
  mn = (n<m ? 0 : n-m);
  mx = (n<l ? n : l);
  ss = 0.0;
  for (i=mn; i<=mx; i++) {
    b1 = math_binom_f(l,i);
    b2 = math_binom_f(m,n-i);
    ss += (b1*b2*pow(a,l-i)*pow(b,m+i-n));
    }
  return(ss);
  }

/* evaluate auxiliary A(l1,l2,A,B,C,g) function
 
   i1,i2,i3 - integer specifiers
   n1,n2    - power IDs for specific Cartesian components of GTOs
   A,B      - specific Cartesian coordinate of the GTOs
   C        - specific Cartesian coordinate of the nucleus
   g        - sum of the GTO exponents */
double gto_afce_a(unsigned i1, unsigned i2, unsigned i3,
  unsigned n1, unsigned n2, double a, double b, double c, double g) {
  static double ff,fu,fd,fa;
  ff = gto_afce_f(i1,n1,n2,a,b);
  fu = math_fact_f(i1)*pow(0.25/g,i2+i3);
  if (i1>=(2*(i2+i3)))
    fu *= pow(c,i1-2*(i2+i3));
  else
    fu /= pow(c,2*(i2+i3)-i1);
  fd = math_fact_f(i2)*math_fact_f(i3);
  if (i1>(2*(i2+i3)))
    fd *= math_fact_f(i1-2*(i2+i3));
  fa = (ff*fu/fd);
  if ((i1+i3)%2)
    fa = (-fa);
  return(fa);
  }

/* -------------------------------------------------------------------------- */
