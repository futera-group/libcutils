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

/* allocate memory for matrix of int-s 

   nrow - number of rows of the matrix
   ncol - number of columns of the matrix */
int **mat_ialloc(unsigned int nrow, unsigned int ncol) {
  int **mat,*p;
  unsigned i;
  mat = (int**)calloc(nrow,sizeof(int*));
  if (!mat)
   msg_error("cannot allocate memory for int matrix",1);
  for (i=0; i<nrow; i++) {
    p = (int*)calloc(ncol,sizeof(int));
    if (!p)
      msg_error("cannot allocate memory for column of int matrix",1);
    mat[i] = p;
    }
  return(mat);
  }

/* allocate memory for matrix of unsigned int-s

   nrow - number of rows of the matrix
   ncol - number of columns of the matrix */
unsigned **mat_ualloc(unsigned nrow, unsigned ncol) {
  unsigned i,**mat,*p;
  mat = (unsigned**)calloc(nrow,sizeof(unsigned*));
  if (!mat)
   msg_error("cannot allocate memory for unsigned int matrix",1);
  for (i=0; i<nrow; i++) {
    p = (unsigned*)calloc(ncol,sizeof(unsigned));
    if (!p)
      msg_error("cannot allocate memory for column of unsigned int matrix",1);
    mat[i] = p;
    }
  return(mat);
  }

/* allocate memory for matrix of short int-s 

   nrow - number of rows of the matrix
   ncol - number of columns of the matrix */
short int **mat_sialloc(unsigned int nrow, unsigned int ncol) {
  short int **mat,*p;
  unsigned i;
  mat = (short int**)calloc(nrow,sizeof(short int*));
  if (!mat)
   msg_error("cannot allocate memory for int matrix",1);
  for (i=0; i<nrow; i++) {
    p = (short int*)calloc(ncol,sizeof(short int));
    if (!p)
      msg_error("cannot allocate memory for column of int matrix",1);
    mat[i] = p;
    }
  return(mat);
  }

/* allocate memory for matrix of unsigned long int-s

   nrow - number of rows of the matrix
   ncol - number of columns of the matrix */
long unsigned **mat_lualloc(unsigned nrow, unsigned ncol) {
  long unsigned i,**mat,*p;
  mat = (long unsigned**)calloc(nrow,sizeof(long unsigned*));
  if (!mat)
   msg_error("cannot allocate memory for unsigned long int matrix",1);
  for (i=0; i<nrow; i++) {
    p = (long unsigned*)calloc(ncol,sizeof(long unsigned));
    if (!p)
      msg_error("cannot allocate memory for column of unsigned long int matrix",1);
    mat[i] = p;
    }
  return(mat);
  }

/* allocate memory for matrix of double precision floating points

   nrow - number of rows of the matrix
   ncol - number of columns of the matrix */
double **mat_falloc(unsigned nrow, unsigned ncol) {
  unsigned i;
  double **mat,*p;
  mat = (double**)calloc(nrow,sizeof(double*));
  if (!mat)
   msg_error("cannot allocate memory for double matrix",1);
  for (i=0; i<nrow; i++) {
    p = (double*)calloc(ncol,sizeof(double));
    if (!p)
      msg_error("cannot allocate memory for column of double matrix",1);
    mat[i] = p;
    }
  return(mat);
  }

/* allocate memory for complex matrix

   nr - number of rows
   nc - number of columns */
double complex** mat_zalloc(unsigned nr, unsigned nc) {
  double complex** x = NULL;
  unsigned i;
  x = (double complex**)malloc(nr*sizeof(double complex*));
  if (!x)
    msg_error("cannot allocate complex matrix",1);
  for (i=0; i<nr; i++) {
    x[i] = NULL;
    x[i] = (double complex*)malloc(nc*sizeof(double complex));
    if (!x[i])
      msg_error("cannot allocate complex matrix",1);
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */

/* free memory of matrix of int-s

   mat  - matrix of doubles
   ncol - number of columns of the matrix */
void* mat_ifree(int **mat, unsigned int ncol) {
  unsigned i;
  if (mat) {
    for (i=0; i<ncol; i++)
      free(mat[i]);
    free(mat);
    }
  return(NULL);
  }

/* free memory of matrix of unsigned int-s

   mat  - matrix of unsigned int-s
   ncol - number of columns of the matrix */
void* mat_ufree(unsigned **mat, unsigned ncol) {
  unsigned i;
  if (mat && ncol) {
    for (i=0; i<ncol; i++)
      free(mat[i]);
    free(mat);
    }
  return(NULL);
  }

/* free memory of matrix of short int-s

   mat  - matrix of doubles
   ncol - number of columns of the matrix */
void* mat_sifree(short int **mat, unsigned int ncol) {
  unsigned i;
  if (mat) {
    for (i=0; i<ncol; i++)
      free(mat[i]);
    free(mat);
    }
  return(NULL);
  }

/* free memory of matrix of unsigned long int-s

   mat  - matrix of unsigned long int-s
   ncol - number of columns of the matrix */
void* mat_lufree(long unsigned **mat, unsigned ncol) {
  unsigned i;
  if (mat && ncol) {
    for (i=0; i<ncol; i++)
      free(mat[i]);
    free(mat);
    }
  return(NULL);
  }

/* free memory of matrix of doubles

   mat  - matrix of doubles
   ncol - number of columns of the matrix */
void* mat_ffree(double **mat, unsigned int ncol) {
  unsigned i;
  if (mat) {
    for (i=0; i<ncol; i++)
      free(mat[i]);
    free(mat);
    }
  return(NULL);
  }

/* free memory allocated for complex matrix

   x  - the matrix
   nr - number of rows */
double complex** mat_zfree(double complex **x, unsigned nr) {
  unsigned i;
  if (x) {
    for (i=0; i<nr; i++)
      if (x[i])
        free(x[i]);
    free(x);
    x = NULL;
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */
