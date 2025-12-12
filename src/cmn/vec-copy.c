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
#include <stdlib.h>
#include <string.h>
#include "cmn/message.h"
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* copy values from one vector to another (integer version)

   vec1 - vector where the values are saved
   vec2 - vector where the values are read from 
   num  - dimension of the vectors */
void vec_icopy(int *vec1, int *vec2, unsigned num) {
  unsigned i;
  if (!vec1 || !vec2)
    msg_error("null pointer in copy int vector",1);
  for (i=0; i<num; i++)
    vec1[i] = vec2[i];
  }

/* create new vector and fill it with values from another one (int version)

   vec - vector where the values are read from 
   num - dimension of the vectors */
int* vec_icopy_new(int *vec, unsigned num) {
  int *v;
  unsigned i;
  if (!vec)
    return(NULL);
  v = vec_ialloc(num);
  for (i=0; i<num; i++)
    v[i] = vec[i];
  return(v);
  }

/* -------------------------------------------------------------------------- */

/* copy values from one vector to another (unsigned integer version)

   vec1 - vector where the values are saved
   vec2 - vector where the values are read from 
   num  - dimension of the vectors */
void vec_ucopy(unsigned *vec1, unsigned *vec2, unsigned num) {
  unsigned i;
  if (!vec1 || !vec2)
    msg_error("null pointer in copy unsigned int vector",1);
  for (i=0; i<num; i++)
    vec1[i] = vec2[i];
  }

/* create new vector and fill it with values from another one (uint version)

   vec - vector where the values are read from 
   num - dimension of the vectors */
unsigned* vec_ucopy_new(unsigned *vec, unsigned num) {
  unsigned *v = NULL;
  unsigned i;
  if (!vec)
    return(NULL);
  v = vec_ualloc(num);
  for (i=0; i<num; i++)
    v[i] = vec[i];
  return(v);
  }

/* -------------------------------------------------------------------------- */

/* copy values from one vector to another (short integer version)

   vec1 - vector where the values are saved
   vec2 - vector where the values are read from 
   num  - dimension of the vectors */
void vec_sicopy(short *vec1, short *vec2, unsigned num) {
  unsigned i;
  if (!vec1 || !vec2)
    msg_error("null pointer in copy s-int vector",1);
  for (i=0; i<num; i++)
    vec1[i] = vec2[i];
  }

/* create new vector and fill it with values from another one (short int)

   vec - vector where the values are read from 
   num - dimension of the vectors */
short* vec_sicopy_new(short *vec, unsigned num) {
  short *v;
  unsigned i;
  if (!vec)
    return(NULL);
  v = vec_sialloc(num);
  for (i=0; i<num; i++)
    v[i] = vec[i];
  return(v);
  }

/* -------------------------------------------------------------------------- */

/* copy values from one vector to another (unsigned long integer version)

   vec1 - vector where the values are saved
   vec2 - vector where the values are read from 
   num  - dimension of the vectors */
void vec_lucopy(unsigned long *vec1, unsigned long *vec2, unsigned num) {
  unsigned i;
  if (!vec1 || !vec2)
    msg_error("null pointer in copy unsigned int vector",1);
  for (i=0; i<num; i++)
    vec1[i] = vec2[i];
  }

/* create new vector and fill it with values from another one (luint version)

   vec - vector where the values are read from 
   num - dimension of the vectors */
unsigned long *vec_lucopy_new(unsigned long *vec, unsigned num) {
  unsigned long *v;
  unsigned i;
  if (!vec)
    return(NULL);
  v = vec_lualloc(num);
  for (i=0; i<num; i++)
    v[i] = vec[i];
  return(v);
  }

/* -------------------------------------------------------------------------- */

/* copy values from one vector to another (double precision real version)

   vec1 - vector where the values are saved
   vec2 - vector where the values are read from 
   num  - dimension of the vectors */
void vec_fcopy(double *vec1, double *vec2, unsigned num) {
  unsigned i;
  if (!vec1 || !vec2)
    msg_error("null pointer in copy double vector",1);
  for (i=0; i<num; i++)
    vec1[i] = vec2[i];
  }

/* copy values from one vector to another (long double precision real version)

   vec1 - vector where the values are saved
   vec2 - vector where the values are read from 
   num  - dimension of the vectors */
void vec_lfcopy(long double *vec1, long double *vec2, unsigned num) {
  unsigned i;
  if (!vec1 || !vec2)
    msg_error("null pointer in copy double vector",1);
  for (i=0; i<num; i++)
    vec1[i] = vec2[i];
  }

/* create new vector and fill it with values from another one (double version)

   vec - vector where the values are read from 
   num - dimension of the vectors */
double* vec_fcopy_new(double *vec, unsigned num) {
  double *v;
  unsigned i;
  if (!vec)
    return(NULL);
  v = vec_falloc(num);
  for (i=0; i<num; i++)
    v[i] = vec[i];
  return(v);
  }

/* copy values from one vector to another and scale them (double)

   vec1 - vector where the values are saved
   vec2 - vector where the values are read from 
   scl  - scaling factor
   num  - dimension of the vectors */
void vec_fcopy_scaled(double *vec1, double *vec2, double scl, unsigned num) {
  unsigned i;
  if (!vec1 || !vec2)
    msg_error("null pointer in copy double vector",1);
  for (i=0; i<num; i++)
    vec1[i] = vec2[i]*scl;
  }

/* -------------------------------------------------------------------------- */

/* copy values from one vector to another (double precision complex version)

   vec1 - vector where the values are saved
   vec2 - vector where the values are read from 
   num  - dimension of the vectors */
void vec_zcopy(double complex *vec1, double complex *vec2, unsigned num) {
  unsigned i;
  if (!vec1 || !vec2)
    msg_error("null pointer in copy complex vector",1);
  for (i=0; i<num; i++)
    vec1[i] = vec2[i];
  }

/* create new vector and fill it with values from another one (complex version)

   vec - vector where the values are read from 
   num - dimension of the vectors */
double complex* vec_zcopy_new(double complex *vec, unsigned num) {
  double complex *v;
  unsigned i;
  if (!vec)
    return(NULL);
  v = vec_zalloc(num);
  for (i=0; i<num; i++)
    v[i] = vec[i];
  return(v);
  }

/* -------------------------------------------------------------------------- */

/* copy values from one vector to another (string version)

   vec1 - vector where the values are saved
   vec2 - vector where the values are read from 
   num  - dimension of the vectors */
void vec_scopy(char **vec1, char **vec2, unsigned num) {
  unsigned i;
  if (!vec1 || !vec2)
    msg_error("null pointer in copy string vector",1);
  for (i=0; i<num; i++)
    strcpy(vec1[i],vec2[i]);
  }

/* create new vector and fill it with values from another one (string version)

   vec   - vector where the values are read from 
   num   - dimension of the vectors
   width - width of single string */
char** vec_scopy_new(char **vec, unsigned num, unsigned width) {
  char **v;
  unsigned i;
  if (!vec)
    return(NULL);
  v = vec_salloc(num,width);
  for (i=0; i<num; i++)
    strcpy(v[i],vec[i]);
  return(v);
  }

/* -------------------------------------------------------------------------- */

/* copy values from one vector to another (user defined type version)

   vec1 - vector where the values are saved
   vec2 - vector where the values are read from 
   f    - function for copying user defined struct
   size - size of the user defined type
   num  - dimension of the vectors */
void vec_tcopy(void *vec1, void *vec2, void(f)(void*,void*), unsigned size,
  unsigned num) {
  unsigned i;
  if (!vec1 || !vec2)
    msg_error("null pointer in copy user defined vector",1);
  for (i=0; i<num; i++) 
    f(vec1+size*i,vec2+size*i);
  }

/* create new vector and fill it with values from another one (user def.)

   vec  - vector where the values are read from 
   f    - function for copying user defined struct
   size - size of one struct in the vector
   num  - dimension of the vectors */
void* vec_tcopy_new(void *vec, void(f)(void*,void*),
  unsigned size, unsigned num) {
  void *v;
  unsigned i;
  if (!vec)
    return(NULL);
  v = vec_talloc(size,num);
  for (i=0; i<num; i++) 
    f(v+size*i,vec+size*i);
  return(v);
  }

/* -------------------------------------------------------------------------- */
