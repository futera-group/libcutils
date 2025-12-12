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
#include "cmn/message.h"
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* change dimension of vector (int)

   v - vector of integers
   n - actual dimension of the vector
   m - new dimension of the vector */
int *vec_iresize(int *v, unsigned n, unsigned m) {
  int *vec;
  unsigned i,min;
  if (n==m)
    return(v);
  if (!v && n)
    msg_error("null pointer as i-vector for resizing",1);
  vec = vec_ialloc(m);
  vec_iset(vec,0,m);
  min = (n<m ? n : m);
  for (i=0; i<min; i++)
    vec[i] = v[i];
  vec_ifree(v);
  return(vec);
  }

/* change dimension of vector (u-int)

   v - vector of unsigned integers
   n - actual dimension of the vector
   m - new dimension of the vector */
unsigned *vec_uresize(unsigned *v, unsigned n, unsigned m) {
  unsigned i,min,*vec;
  if (n==m)
    return(v);
  if (!v && n)
    msg_error("null pointer as u-vector for resizing",1);
  vec = vec_ualloc(m);
  vec_uset(vec,0,m);
  min = (n<m ? n : m);
  for (i=0; i<min; i++)
    vec[i] = v[i];
  vec_ufree(v);
  return(vec);
  }

/* change dimension of vector (s-int)

   v - vector of short integers
   n - actual dimension of the vector
   m - new dimension of the vector */
short *vec_siresize(short *v, unsigned n, unsigned m) {
  short *vec;
  unsigned i,min;
  if (n==m)
    return(v);
  if (!v && n)
    msg_error("null pointer as si-vector for resizing",1);
  vec = vec_sialloc(m);
  vec_siset(vec,0,m);
  min = (n<m ? n : m);
  for (i=0; i<min; i++)
    vec[i] = v[i];
  vec_sifree(v);
  return(vec);
  }

/* change dimension of vector (l-int)

   v - vector of long integers
   n - actual dimension of the vector
   m - new dimension of the vector */
long *vec_liresize(long *v, unsigned n, unsigned m) {
  unsigned i,min;
  long *vec;
  if (n==m)
    return(v);
  if (!v && n)
    msg_error("null pointer as li-vector for resizing",1);
  vec = vec_lialloc(m);
  vec_liset(vec,0,m);
  min = (n<m ? n : m);
  for (i=0; i<min; i++)
    vec[i] = v[i];
  vec_lifree(v);
  return(vec);
  }

/* change dimension of vector (lu-int)

   v - vector of long unsigned integers
   n - actual dimension of the vector
   m - new dimension of the vector */
long unsigned *vec_luresize(long unsigned*v, unsigned n, unsigned m) {
  long unsigned *vec;
  unsigned i,min;
  if (n==m)
    return(v);
  if (!v && n)
    msg_error("null pointer as lu-vector for resizing",1);
  vec = vec_lualloc(m);
  vec_luset(vec,0,m);
  min = (n<m ? n : m);
  for (i=0; i<min; i++)
    vec[i] = v[i];
  vec_lufree(v);
  return(vec);
  }

/* change dimension of vector (double)

   v - vector of double precision reals
   n - actual dimension of the vector
   m - new dimension of the vector */
double *vec_fresize(double *v, unsigned n, unsigned m) {
  double *vec;
  unsigned i,min;
  if (n==m)
    return(v);
  if (!v && n)
    msg_error("null pointer as f-vector for resizing",1);
  vec = vec_falloc(m);
  vec_fset(vec,0.0,m);
  min = (n<m ? n : m);
  for (i=0; i<min; i++)
    vec[i] = v[i];
  vec_ffree(v);
  return(vec);
  }

/* change dimension of vector (string)

   v - vector of integers
   n - actual dimension of the vector
   m - new dimension of the vector */
char **vec_sresize(char **s, unsigned n, unsigned m) {
  char **vec;
  unsigned i;
  if (n==m)
    return(s);
  if (!s && n)
    msg_error("null pointer as s-vector for resizing",1);
  vec = vec_talloc(sizeof(char*),m);
  for (i=0; i<m; i++) {
    vec[i] = NULL;
    if (i<n) {
      vec[i] = s[i];
      s[i] = NULL;
      }
    }
  vec_sfree(s,n);
  return(vec);
  }

/* change dimension of vector (user defined type)

   v - vector of user defined types
   n - actual dimension of the vector
   m - new dimension of the vector
   s - size of one element of the vector
   f - function for copying the structs */
void *vec_tresize(void *v, unsigned n, unsigned m, unsigned s, 
  void(f)(void*,void*)) {
  unsigned i,min;
  void *vec;
  if (n==m)
    return(v);
  if (!v && n) 
    msg_error("null pointer as t-vector for resizing",1);
  vec = vec_talloc(s,m);
  min = (n<m ? n : m);
  for (i=0; i<min; i++)
    f(vec+s*i,v+s*i);
  vec_tfree(v);
  return(vec);
  }

/* -------------------------------------------------------------------------- */
