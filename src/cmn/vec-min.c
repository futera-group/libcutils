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
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* find the smallest value in given vector (unsinged int version)

   v - unsigned int vector
   n - number of elements in the vector */
unsigned vec_umin(unsigned *v, unsigned n) {
  unsigned i,id,min;
  if (n==0)
    return(0);
  id = 0;
  min = v[0];
  for (i=0; i<n; i++)
    if (min>v[i]) {
      id = i;
      min = v[i];
      }
  return(id);
  }

/* find the largest value in given vector (unsinged int version)

   v - unsigned int vector
   n - number of elements in the vector */
unsigned vec_umax(unsigned *v, unsigned n) {
  unsigned i,id,max;
  if (n==0)
    return(0);
  id = 0;
  max = v[0];
  for (i=0; i<n; i++)
    if (max<v[i]) {
      id = i;
      max = v[i];
      }
  return(id);
  }

/* -------------------------------------------------------------------------- */

/* find minimum value in given vector (double precision real version)

   v - double precision vector
   n - number of elements in the vector */
unsigned vec_fmin(double *v, unsigned n) {
  unsigned i,id;
  double min;
  if (n==0)
    return(0);
  id = 0;
  min = v[0];
  for (i=0; i<n; i++)
    if (min>v[i]) {
      id = i;
      min = v[i];
      }
  return(id);
  }

/* find minimum abs value in given vector (double precision real version) 

   v - double precision vector
   n - number of elements in the vector */
unsigned vec_fmin_abs(double *v, unsigned n) {
  unsigned i,id;
  double min,val;
  if (n==0)
    return(0);
  id = 0;
  min = 0.0;
  for (i=0; i<n; i++) {
    val=fabs(v[i]);
    if (min>val) {
      id = i;
      min = val;
      }
    }
  return(id);
  }

/* find minimum and maximum from values in the vector (double)
 
   v  - double precision vector
   n  - number of values in the vector
   mn - lowest value
   mx - highest value */
void vec_fminmax(double *v, long unsigned n, double *mn, double *mx) {
  double min,max;
  long unsigned i;
  if (n) {
    min = max = v[0];
    for (i=1; i<n; i++) {
      if (v[i]<min)
        min = v[i];
      if (v[i]>max)
        max = v[i];
      }
    (*mn) = min;
    (*mx) = max;
    }
  else {
    (*mn) = 0.0;
    (*mx) = 0.0;
    }
  }

/* find maximum value in given vector (double precision real version)

   v - double precision vector
   n - number of elements in the vector */
unsigned vec_fmax(double *v, unsigned n) {
  unsigned i,id;
  double max;
  if (n==0)
    return(0);
  id = 0;
  max = v[0];
  for (i=0; i<n; i++)
    if (max<v[i]) {
      id = i;
      max = v[i];
      }
  return(id);
  }

/* find maximum abs value in given vector (double precision real version)

   v - double precision vector
   n - number of elements in the vector */
unsigned vec_fmax_abs(double *v, unsigned n) {
  unsigned i,id;
  double max,val;
  if (n==0)
    return(0);
  id = 0;
  max = 0.0;
  for (i=0; i<n; i++) {
    val = fabs(v[i]);
    if (max<val) {
      id = i;
      max = val;
      }
    }
  return(id);
  }

/* -------------------------------------------------------------------------- */

/* find longest string  in given vector (string version)
 
   vec - array of strings
   num - length of the array */
unsigned vec_smax(char **vec, unsigned num) {
  unsigned i;
  size_t n,max;
  max = 0;
  for (i=0; i<num; i++) {
    n = str_length(vec[i]);
    if (n>max)
      max = n;
    }
  return(max);
  }

/* -------------------------------------------------------------------------- */
