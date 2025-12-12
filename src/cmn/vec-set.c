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

#include <complex.h>
#include <string.h>

/* -------------------------------------------------------------------------- */

/* set the whole integer vector on one value

   vec - vector of int-s
   n   - dimension of the vector
   val - the value for setting */
void vec_iset(int *vec, int val, long unsigned n) {
  long unsigned i;
  for (i=0; i<n; i++)
    vec[i] = val;
  }

/* set the whole unsigned integer vector on one value

   vec - vector of unsigned int-s
   n   - dimension of the vector
   val - the value for setting */
void vec_uset(unsigned *vec, unsigned val, long unsigned n) {
  long unsigned i;
  for (i=0; i<n; i++)
    vec[i] = val;
  }

/* set the whole short integer vector on one value

   vec - vector of short int-s
   n   - dimension of the vector
   val - the value for setting */
void vec_siset(short *vec, short val, long unsigned n) {
  long unsigned i;
  for (i=0; i<n; i++)
    vec[i] = val;
  }

/* set the whole long integer vector on one value

   vec - vector of long int-s
   n   - dimension of the vector
   val - the value for setting */
void vec_liset(long *vec, long val, long unsigned n) {
  long unsigned i;
  for (i=0; i<n; i++)
    vec[i] = val;
  }

/* set the whole unsigned long integer vector on one value

   vec - vector of unsigned long int-s
   n   - dimension of the vector
   val - the value for setting */
void vec_luset(unsigned long *vec, unsigned long val, long unsigned n) {
  long unsigned i;
  for (i=0; i<n; i++)
    vec[i] = val;
  }

/* set the whole floating point vector on one value

   vec - vector of reals
   n   - dimension of the vector
   val - the value for setting */
void vec_fset(double *vec, double val, long unsigned n) {
  long unsigned i;
  for (i=0; i<n; i++)
    vec[i] = val;
  }

/* set the whole floating point vector on one value

   vec - vector of reals
   n   - dimension of the vector
   val - the value for setting */
void vec_lfset(long double *vec, long double val, long unsigned n) {
  long unsigned i;
  for (i=0; i<n; i++)
    vec[i] = val;
  }

/* set the whole complex real vector on one value

   vec - vector of reals
   n   - dimension of the vector
   val - the value for setting */
void vec_zset(complex double *vec, complex double val, long unsigned n) {
  long unsigned i;
  for (i=0; i<n; i++)
    vec[i] = val;
  }

/* set the whole string vector on one value

   vec - vector of strings
   n   - dimension of the vector
   val - the value for setting */
void vec_sset(char **vec, char *val, long unsigned n) {
  long unsigned i;
  for (i=0; i<n; i++)
    strcpy(vec[i],val);
  }

/* -------------------------------------------------------------------------- */
