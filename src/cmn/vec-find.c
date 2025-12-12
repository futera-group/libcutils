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

#include <string.h>
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* search given value in vector (int version)

   vec - vector of int-s
   val - the value which is searched
   num - dimension of the vector */
int vec_ifind(int *vec, int val, unsigned num) {
  unsigned i;
  if (!vec)
    return(-1);
  for (i=0; i<num; i++)
    if (vec[i]==val)
      return(i);
  return(-1);
  }

/* search given value in vector (unsigned int version)

   vec - vector of unsigned int-s
   val - the value which is searched
   num - dimension of the vector */
int vec_ufind(unsigned *vec, unsigned val, unsigned num) {
  unsigned i;
  if (!vec)
    return(-1);
  for (i=0; i<num; i++)
    if (vec[i]==val)
      return(i);
  return(-1);
  }

/* search given value in vector (short int version)

   vec - vector of short int-s
   val - the value which is searched
   num - dimension of the vector */
int vec_sifind(short *vec, short val, unsigned num) {
  unsigned i;
  if (!vec)
    return(-1);
  for (i=0; i<num; i++)
    if (vec[i]==val)
      return(i);
  return(-1);
  }

/* search given value in vector (unsigned long int version)

   vec - vector of unsigned long int-s
   val - the value which is searched
   num - dimension of the vector */
int vec_lufind(unsigned long *vec, unsigned long val, unsigned num) {
  unsigned i;
  if (!vec)
    return(-1);
  for (i=0; i<num; i++)
    if (vec[i]==val)
      return(i);
  return(-1);
  }

/* search given value in vector (real version)

   vec - vector of double-s
   val - the value which is searched
   num - dimension of the vector */
int vec_ffind(double *vec, double val, unsigned num) {
  unsigned i;
  if (!vec)
    return(-1);
  for (i=0; i<num; i++)
    if (vec[i]==val)
      return(i);
  return(-1);
  }

/* search given value in vector (string version)

   vec - vector of strings
   num - dimension of the vector
   val - the value which is searched */
int vec_sfind(char **vec, char *val, unsigned num) {
  unsigned i;
  if (!vec)
    return(-1);
  for (i=0; i<num; i++)
    if (str_compare(vec[i],val))
      return(i);
  return(-1);
  }

/* search given value in vector (user defined types)

   vec  - vector of user defined types
   val  - the value which is searched
   f    - function for comparing user defined types
   size - size of the user defined type
   num  - dimension of the vector */
int vec_tfind(void *vec, void *val, int(f)(void*,void*), unsigned size,
  unsigned num) {
  unsigned i;
  if (!vec)
    return(-1);
  for (i=0; i<num; i++)
    if (f(vec+i*size,val))
      return(i);
  return(-1);
  }

/* -------------------------------------------------------------------------- */
