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
#include <cmn/message.h>
#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* calculate normalization constant for Cartesian GTO

   e - exponent of the function
   n - GTO specifier (l,n,m) */
double gto_norm(double e, unsigned *n) {
  static double t=1.0,b=1.0;
  static unsigned a=1,s;
  s = n[0]+n[1]+n[2];
  switch (s) {
    case 0: /* S fce */
      t = 8.0;
      b = 1.0;
      a = 3;
      break;
    case 1: /* P fce */
      t = 128.0;
      b = 1.0;
      a = 5;
      break;
    case 2: /* D fce */
      t = 2048.0;
      a = 7;
      b=((n[0]==2 || n[1]==2 || n[2]==2) ? 9.0 : 1.0);
      break;
    case 3: /* F fce */
      t = 32768.0;
      a = 9;
      if (n[0]==1 && n[1]==1 && n[2]==1)
        b = 1.0;
      else if (n[0]==3 || n[1]==3 || n[2]==3)
        b = 225.0;
      else
        b = 9.0;
      break;
    case 4: /* G fce */
      t = 524288.0; 
      a = 11;
      if (n[0]==4 || n[1]==4 || n[2]==4)
        b = 11025.0;
      else if (n[0]==3 || n[1]==3 || n[2]==3)
        b = 225.0;
      else if (n[0]==1 || n[1]==1 || n[2]==1)
        b = 9.0;
      else
        b = 81.0;
      break;
    case 5: /* H fce */
      t = 8388608.0;
      a = 13;
      if (n[0]==5 || n[1]==5 || n[2]==5)
        b = 893025.0;
      else if (n[0]==4 || n[1]==4 || n[2]==4)
        b = 11025.0;
      else if (n[0]==3 || n[1]==3 || n[2]==3) {
        if (n[0]==2 || n[1]==2 || n[3]==2)
          b = 2025.0;
        else
          b = 225.0;
        }
      else
        b = 81.0;
      break;
    default:
      msg_error_f("unsupported type in gauss_basis_bfce_norm (s = %d)",1,s);
    }
  /* return normalization factor */
  return(pow((t*pow(e,a))/(b*MATH_PIs3),0.25));
  }

/* -------------------------------------------------------------------------- */
