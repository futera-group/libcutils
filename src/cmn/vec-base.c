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
#include "cmn/message.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for vector of int-s
   
   amax - length of the vector */
int *vec_ialloc(long unsigned vmax) {
  int *vec = NULL;
  if (vmax) {
    vec = calloc(vmax,sizeof(int));
    if (!vec)
      msg_error("cannot allocate memory for vector of int-s",1);
    }
  return(vec);
  }

/* allocate memory for vector of unsigned int-s
   
   amax - length of the vector */
unsigned *vec_ualloc(long unsigned vmax) {
  unsigned *vec=NULL;
  if (vmax) {
    vec = calloc(vmax,sizeof(unsigned));
    if (!vec)
      msg_error("cannot allocate memory for vector of unsigned int-s",1);
    }
  return(vec);
  }

/* allocate memory for vector of short int-s
   
   amax - length of the vector */
short *vec_sialloc(long unsigned vmax) {
  short *vec = NULL;
  if (vmax) {
    vec = calloc(vmax,sizeof(short));
    if (!vec)
      msg_error("cannot allocate memory for vector of short int-s",1);
    }
  return(vec);
  }

/* allocate memory for vector of short unsigned int-s
   
   amax - length of the vector */
short unsigned *vec_sualloc(long unsigned vmax) {
  short unsigned *vec = NULL;
  if (vmax) {
    vec = calloc(vmax,sizeof(short unsigned));
    if (!vec)
      msg_error("cannot allocate memory for vector of short unsigned int-s",1);
    }
  return(vec);
  }

/* allocate memory for vector of long integers
   
   amax - length of the vector */
long *vec_lialloc(long unsigned vmax) {
  long *vec = NULL;
  if (vmax) {
    vec = calloc(vmax,sizeof(long));
    if (!vec)
      msg_error("cannot allocate memory for vector of long integers",1);
    }
  return(vec);
  }

/* allocate memory for vector of unsigned longs
   
   amax - length of the vector */
unsigned long *vec_lualloc(long unsigned vmax) {
  unsigned long *vec = NULL;
  if (vmax) {
    vec = calloc(vmax,sizeof(unsigned long));
    if (!vec)
      msg_error("cannot allocate memory for vector of unsigned longs",1);
    }
  return(vec);
  }

/* allocate memory for array of strings
   
   amax - length of the array
   smax - length of strings */
char **vec_salloc(long unsigned amax, long unsigned smax) {
  long unsigned i;
  char **vec = NULL;
  if (amax && smax) {
    vec = (char**)malloc(amax*sizeof(char*));
    if (!vec)
      msg_error("cannot allocate memory for array of strings",1);
    for (i=0; i<amax; i++) {
      vec[i] = NULL;
      vec[i] = (char*)malloc(smax*sizeof(char));
      if (!vec[i])
        msg_error("cannot allocate memory for string in array",1);
      }
    }
  return(vec);
  }

/* allocate memory for vector of double precision reals */
double *vec_falloc(long unsigned vmax) {
  double *vec = NULL;
  if (vmax) {
    vec = calloc(vmax,sizeof(double));
    if (!vec)
      msg_error("cannot allocate memory for vector of doubles",1);
    }
  return(vec);
  }

/* allocate memory for vector of longd double precision reals */
long double *vec_lfalloc(long unsigned vmax) {
  long double *vec = NULL;
  if (vmax) {
    vec = calloc(vmax,sizeof(long double));
    if (!vec)
      msg_error("cannot allocate memory for vector of long doubles",1);
    }
  return(vec);
  }

/* allocate memory for vector of double precision complex numbers */
double complex *vec_zalloc(long unsigned vmax) {
  double complex *vec = NULL;
  if (vmax) {
    vec = calloc(vmax,sizeof(double complex));
    if (!vec)
      msg_error("cannot allocate memory for vector of complex numbers",1);
    }
  return(vec);
  }

/* allocate memory for array of user defined structs
   
   asize - size of one item in the array
   amax  - length of the array */
void *vec_talloc(long unsigned asize, long unsigned amax) {
  void *vec = NULL;
  if (asize && amax) {
    vec = malloc(asize*amax);
    if (!vec)
      msg_error("cannot allocate memory for user defined vector",1);
    }
  return(vec);
  }

/* -------------------------------------------------------------------------- */

/* deallocate memory of vector of int-s
 
   vec - vector of int-s */
void* vec_ifree(int *vec) {
  if (vec)
    free(vec);
  return(NULL);
  }

/* deallocate memory of vector of unsigned int-s
 
   vec - vector of unsigned int-s */
void* vec_ufree(unsigned *vec) {
  if (vec)
    free(vec);
  return(NULL);
  }

/* deallocate memory of vector of short int-s
 
   vec - vector of short int-s */
void* vec_sifree(short *vec) {
  if (vec)
    free(vec);
  return(NULL);
  }

/* deallocate memory of vector of short unsigned int-s
 
   vec - vector of short int-s */
void* vec_sufree(short unsigned *vec) {
  if (vec)
    free(vec);
  return(NULL);
  }

/* deallocate memory of vector of unsigned longs
 
   vec - vector of unsigned longs */
void* vec_lifree(long *vec) {
  if (vec)
    free(vec);
  return(NULL);
  }

/* deallocate memory of vector of unsigned longs
 
   vec - vector of unsigned longs */
void* vec_lufree(unsigned long *vec) {
  if (vec)
    free(vec);
  return(NULL);
  }

/* deallocate memory of vector of double precision reals */
void* vec_ffree(double *vec) {
  if (vec)
    free(vec);
  return(NULL);
  }

/* deallocate memory of vector of long double precision reals */
void* vec_lffree(long double *vec) {
  if (vec)
    free(vec);
  return(NULL);
  }

/* deallocate memory of vector of double precision complex numbers */
void* vec_zfree(double complex *vec) {
  if (vec)
    free(vec);
  return(NULL);
  }

/* deallocate memory of array of strings
 
   array - array of strings
   amax - length of the array */
void* vec_sfree(char **array, long unsigned amax) {
  long unsigned i;
  if (!array)
    return(NULL);
  for (i=0; i<amax; i++)
    free(array[i]);
  free(array);
  return(NULL);
  }

/* deallocate memory of array of user defined structs
 
   array - array of structs */
void* vec_tfree(void *array) {
  if (array)
    free(array);
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
