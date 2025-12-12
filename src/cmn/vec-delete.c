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

#include <stdlib.h>
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* delete one element in vector (integer)

   v - vector of integers
   n - actual dimension of the vector
   d - ID of the element for deletion */
int *vec_idelete(int *v, unsigned n, unsigned d) {
  int *w=NULL;
  unsigned i;
  if (n>1) {
    w = vec_ialloc(n-1);
    for (i=0; i<d; i++)
      w[i] = v[i];
    for (i=d; i<n-1; i++)
      w[i] = v[i+1];
    vec_ifree(v);
    }
  else
    vec_ifree(v);
  return(w);
  }

/* delete one element in vector (unsigned integer)

   v - vector of double precision reals
   n - actual dimension of the vector
   d - ID of the element for deletion */
unsigned *vec_udelete(unsigned *v, unsigned n, unsigned d) {
  unsigned i,*w = NULL;
  if (n>1) {
    w = vec_ualloc(n-1);
    for (i=0; i<d; i++)
      w[i] = v[i];
    for (i=d; i<n-1; i++)
      w[i] = v[i+1];
    vec_ufree(v);
    }
  else
    vec_ufree(v);
  return(w);
  }

/* delete one element in vector (double)

   v - vector of double precision reals
   n - actual dimension of the vector
   d - ID of the element for deletion */
double *vec_fdelete(double *v, unsigned n, unsigned d) {
  double *w = NULL;
  unsigned i;
  if (n>1) {
    w = vec_falloc(n-1);
    for (i=0; i<d; i++)
      w[i] = v[i];
    for (i=d; i<n-1; i++)
      w[i] = v[i+1];
    vec_ffree(v);
    }
  else
    vec_ffree(v);
  return(w);
  }

/* delete one element in vector (user defined type)

   v - vector of user defined types
   n - actual dimension of the vector
   d - ID of the element for deletion
   s - size of one element in the vector */
void *vec_tdelete(void *v, unsigned n, unsigned d, unsigned s,
  void(f)(void*,void*)) {
  unsigned i;
  void *w = NULL;
  if (n>1) {
    w = vec_talloc(s,n-1);
    for (i=0; i<d; i++)
      f(w+s*i,v+s*i);
    for (i=d; i<n-1; i++)
      f(w+s*i,v+s*(i+1));
    vec_tfree(v);
    }
  else
    vec_tfree(v);
  return(w);
  }

/* -------------------------------------------------------------------------- */
