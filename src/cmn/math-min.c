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

/* -------------------------------------------------------------------------- */

/* return greater of the two values (int) */
int math_imax(int a, int b) {
  return(a>b ? a : b);
  }

/* return smaller of the two values (int)*/
int math_imin(int a, int b) {
  return(a<b ? a : b);
  }

/* return greater of the two values (unsigned) */
unsigned math_umax(unsigned a, unsigned b) {
  return(a>b ? a : b);
  }

/* return smaller of the two values (unsigned)*/
unsigned math_umin(unsigned a, unsigned b) {
  return(a<b ? a : b);
  }

/* return greater of the two values (double) */
double math_fmax(double a, double b) {
  return(a>b ? a : b);
  }

/* return smaller of the two values (double) */
double math_fmin(double a, double b) {
  return(a<b ? a : b);
  }

/* -------------------------------------------------------------------------- */
