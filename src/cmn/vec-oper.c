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
#include <stdlib.h>
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

#define PRECISION 1.0E-10

/* -------------------------------------------------------------------------- */

/* return norm of the floating point vector

   vec - vector of double precision reals
   n   - dimension of the vector */
double vec_fnorm(double *vec, unsigned n) {
  unsigned i;
  double res = 0.0;
  for (i=0; i<n; i++)
    res += vec[i]*vec[i];
  return(sqrt(res));
  }

/* -------------------------------------------------------------------------- */

/* normalize vector to unit (double precision real version) 

   vec - pointer to double real vector 
   num - dimension of the vector */
void vec_funit(double *vec, unsigned num) {
  unsigned i;
  double norm = 0.0;
  for (i=0; i<num; i++)
    norm += vec[i]*vec[i];
  norm = sqrt(norm);
  if (norm>PRECISION)
    for (i=0; i<num; i++)
      vec[i] /= norm;
  }

/* create new normalized vector form old one (double precision real version) 

   vec1 - result unit vector
   vec2 - vector which is normalized
   num  - dimension of the vector */
void vec_funit_copy(double *vec1, double *vec2, unsigned num) {
  unsigned i;
  double norm = 0.0;
  for (i=0; i<num; i++)
    norm += vec2[i]*vec2[i];
  norm = sqrt(norm);
  if (norm>PRECISION) 
    for (i=0; i<num; i++)
      vec1[i] = vec2[i]/norm;
  }

/* -------------------------------------------------------------------------- */

/* scalar product of two floating point vectors
   
   vec1,vec2 - double precision real vectors
   n         - dimension of the vectors */
double vec_fprod_scalar(double *vec1, double *vec2, unsigned n) {
  unsigned i;
  double res=0.0;
  for (i=0; i<n; i++)
    res += vec1[i]*vec2[i];
  return(res);
  }

/* vector product of two floating point vectors
   
   vec1,vec2 - double precision real vectors
   res       - final vector */
void vec_fprod_vector(double *vec1, double *vec2, double *res) {
  res[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
  res[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
  res[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];
  }

/* mixed product of three floating point vectors

   vec1,vec2,vec3 - double precision real vectors */
double vec_fprod_mixed(double *vec1, double *vec2, double *vec3) {
  double v[3];
  vec_fprod_vector(vec2,vec3,v);
  return(vec_fprod_scalar(vec1,v,3));
  }

/* -------------------------------------------------------------------------- */

/* angle between two floating point vectors
   
   vec1,vec2 - double precision real vectors
   n         - dimension of the vectors */
double vec_fangle(double *vec1, double *vec2, unsigned n) {
  double res;
  res = vec_fprod_scalar(vec1,vec2,n);
  res = res/(vec_fnorm(vec1,n)*vec_fnorm(vec2,n));
  if (res<-1.0 && fabs(res+1.0)<PRECISION)
    res = -1.0;
  if (res>1.0 && fabs(res-1.0)<PRECISION)
    res = 1.0;
  return(acos(res));
  }

/* -------------------------------------------------------------------------- */

/* sum two vectors (double precision real version)

   v     - resulting vector
   v1,v2 - vectors which are summed
   n     - dimension of the vectors */
void vec_fadd(double *v, double *v1, double *v2, unsigned n) {
  unsigned i;
  for (i=0; i<n; i++)
    v[i] = v1[i]+v2[i];
  }

/* scale and add one vector to another (double precision real version)

   vec1 - the result vector
   vec2 - vector which is scaled
   fac  - scaling factor
   num  - dimension of the vectors */
void vec_fadd_scaled(double *vec1, double *vec2, double fac, unsigned num) {
  unsigned i;
  for (i=0; i<num; i++)
    vec1[i] += fac*vec2[i];
  }

/* subtract two vectors (double precision real version)

   v     - resulting vector
   v1,v2 - vectors which are summed
   n     - dimension of the vectors */
void vec_fsub(double *v, double *v1, double *v2, unsigned n) {
  unsigned i;
  for (i=0; i<n; i++)
    v[i] = v1[i]-v2[i];
  }

/* return sum of all elements in the vector (double precision real version)
 
   v - vector of double precision reals
   n - number of elements in the vector */
double vec_fsum(double *v, unsigned n) {
  unsigned i;
  double res=0.0;
  for (i=0; i<n; i++)
    res += v[i];
  return(res);
  }

/* return sum of absolute values of all elements (double real version)
 
   v - vector of double precision reals
   n - number of elements in the vector */
double vec_fsum_abs(double *v, unsigned n) {
  unsigned i;
  double res=0.0;
  for (i=0; i<n; i++)
    res += fabs(v[i]);
  return(res);
  }

/* scale all vector elements by given factor (double precision real version)

   vec - vector of doubles
   fac - scaling factor
   num - dimension of the vector */
void vec_fscale(double *vec, double fac, unsigned num) {
  unsigned i;
  for (i=0; i<num; i++)
    vec[i] *= fac;
  }

/* -------------------------------------------------------------------------- */

/* merge two unsigned integer vector to the first one

   v1,v2 - pointers to the vectors
   n1,n2 - number of elements in the vectors */
void vec_umerge(unsigned **v1, unsigned *n1, unsigned *v2, unsigned n2) {
  unsigned i,*w=NULL;
  if (!n2)
    return;
  else if (!n1) {
    (*v1) = vec_ucopy_new(v2,n2);
    (*n1) = n2;
    }
  else {
    w = vec_ualloc(*n1+n2);
    for (i=0; i<(*n1); i++)
      w[i] = (*v1)[i];
    for (i=0; i<n2; i++)
      w[i+(*n1)] = v2[i];
    vec_ufree(*v1);
    (*v1) = w;
    (*n1) += n2;
    }
  }

/* merge two unsigned integer vector to new one

   v1,v2 - pointers to the vectors
   n1,n2 - number of elements in the vectors */
unsigned *vec_umerge_new(unsigned *v1, unsigned n1, unsigned *v2, unsigned n2) {
  unsigned i,*w=NULL;
  if (!n1)
    w = v2;
  else if (!n2)
    w = v1;
  else {
    w = vec_ualloc(n1+n2);
    for (i=0; i<n1; i++)
      w[i] = v1[i];
    for (i=0; i<n2; i++)
      w[i+n1] = v2[i];
    }
  return(w);
  }

/* -------------------------------------------------------------------------- */
