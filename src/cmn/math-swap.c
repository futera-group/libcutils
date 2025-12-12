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

/* swap values between two variables (integers) */
void math_iswap(int *a, int *b) {
  int c;
  c = (*a);
  (*a) = (*b);
  (*b) = c;
  }

/* swap values between two variables (unsigned int-s) */
void math_uswap(unsigned *a, unsigned *b) {
  unsigned c;
  c = (*a);
  (*a) = (*b);
  (*b) = c;
  }

/* swap values between two variables (short int-s) */
void math_siswap(short *a, short *b) {
  short c;
  c = (*a);
  (*a) = (*b);
  (*b) = c;
  }

/* swap values between two variables (doubles) */
void math_fswap(double *a, double *b) {
  double c;
  c = (*a);
  (*a) = (*b);
  (*b) = c;
  }

/* -------------------------------------------------------------------------- */
