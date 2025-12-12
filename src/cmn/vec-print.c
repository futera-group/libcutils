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

#include <stdio.h>
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* raw print of integer vector

   file - pointer to open file
   vec  - vector of integers
   num  - number of vector elements */
void vec_iprint(FILE *file, int *vec, unsigned num) {
  unsigned i;
  if (num) {
    fprintf(file,"%d",vec[0]);
    for (i=1; i<num; i++)
      fprintf(file," %d",vec[i]);
    fprintf(file,"\n");
    }
  }

/* formatted printing of integer vector

   file  - pointer to open file
   vec   - vector of int-s
   ntot  - dimension of the vector
   nrow  - number of int-s in one line
   width - space for one number */
void vec_ifprint(FILE *file, int *vec, unsigned ntot,
  unsigned nrow, unsigned width) {
  unsigned i;
  for (i=0; i<ntot; i++) {
    fprintf(file,"%*d",width,vec[i]);
    if ((i+1)%nrow==0)
      fprintf(file,"\n");
    }
  if (ntot%nrow)
    fprintf(file,"\n");
  }

/* -------------------------------------------------------------------------- */

/* raw print of short-integer vector

   file - pointer to open file
   vec  - vector of short integers
   num  - number of vector elements */
void vec_siprint(FILE *file, short *vec, unsigned num) {
  unsigned i;
  if (num) {
    fprintf(file,"%d",vec[0]);
    for (i=1; i<num; i++)
      fprintf(file," %d",vec[i]);
    fprintf(file,"\n");
    }
  }

/* formatted printing of short-integer vector

   file  - pointer to open file
   vec   - vector of short integers
   ntot  - dimension of the vector
   nrow  - number of short integers in one line
   width - space for one number */
void vec_sifprint(FILE *file, short *vec, unsigned ntot,
  unsigned nrow, unsigned width) {
  unsigned i;
  for (i=0; i<ntot; i++) {
    fprintf(file,"%*d",width,vec[i]);
    if ((i+1)%nrow==0)
      fprintf(file,"\n");
    }
  if (ntot%nrow)
    fprintf(file,"\n");
  }

/* -------------------------------------------------------------------------- */

/* raw print of unsigned integer vector

   file - pointer to open file
   vec  - vector of unsigned integers
   num  - number of vector elements */
void vec_uprint(FILE *file, unsigned *vec, unsigned num) {
  unsigned i;
  if (num) {
    fprintf(file,"%d",vec[0]);
    for (i=1; i<num; i++)
      fprintf(file," %d",vec[i]);
    fprintf(file,"\n");
    }
  }

/* formatted printing of unsigned integer vector

   file  - pointer to open file
   vec   - vector of unsigned int-s
   ntot  - dimension of the vector
   nrow  - number of unsigned int-s in one line
   width - space for one number */
void vec_ufprint(FILE *file, unsigned *vec, unsigned ntot,
  unsigned nrow, unsigned width) {
  unsigned i;
  for (i=0; i<ntot; i++) {
    fprintf(file,"%*d",width,vec[i]);
    if ((i+1)%nrow==0)
      fprintf(file,"\n");
    }
  if (ntot%nrow)
    fprintf(file,"\n");
  }

/* -------------------------------------------------------------------------- */

/* print real vector

   v - real vector
   n - dimension of the vector */
void vec_fprint(double *v, unsigned n) {
  unsigned i;
  for (i=0; i<n; i++)
    printf("%15e",v[i]);
  printf("\n");
  }

/* print double-precision real vector with specific number format
    
   v - the vector
   n - dimension of the vector
   f - the number format */
void vec_fprint_f(double *v, unsigned n, char *f) {
  unsigned i;
  for (i=0; i<n; i++)
    printf(f,v[i]);
  printf("\n");
  }

/* print scaled double precision real vector

   v - real vector
   s - scaling factor
   n - dimension of the vector */
void vec_fprint_scale(double *v, double s, unsigned n) {
  unsigned i;
  for (i=0; i<n; i++)
    printf("%15e",v[i]*s);
  printf("\n");
  }

/* print double precision vector

   v       - double precision real vector
   id1,id2 - first and last element for printing */
void vec_fprint_sub(double *v, unsigned id1, unsigned id2) {
  unsigned i;
  for (i=id1; i<id2; i++)
    printf("%15e",v[i]);
  printf("\n");
  }

/* formatted printing of double precision float vector

   file  - pointer to open file
   vec   - vector of int-s
   ntot  - dimension of the vector
   nrow  - number of int-s in one line
   width - space for one number
   prec  - precision for number printing 
   frm   - format of float number (f,F,e,E,g,G) */
void vec_ffprint(FILE *file, double *vec, unsigned ntot,
  unsigned nrow, unsigned width, unsigned prec, short frm) {
  unsigned i;
  for (i=0; i<ntot; i++) {
    switch (frm) {
      case VEC_FRM_f: fprintf(file,"%*.*f",width,prec,vec[i]); break;
      case VEC_FRM_F: fprintf(file,"%*.*F",width,prec,vec[i]); break;
      case VEC_FRM_e: fprintf(file,"%*.*e",width,prec,vec[i]); break;
      case VEC_FRM_E: fprintf(file,"%*.*E",width,prec,vec[i]); break;
      case VEC_FRM_g: fprintf(file,"%*.*g",width,prec,vec[i]); break;
      case VEC_FRM_G: fprintf(file,"%*.*G",width,prec,vec[i]); break;
      default: fprintf(file,"%*.*f",width,prec,vec[i]); break;
      }
    if ((i+1)%nrow==0)
      fprintf(file,"\n");
    }
  if (ntot%nrow)
    fprintf(file,"\n");
  }

/* formatted printing of double precision float vector to table format

   f  - pointer to open file
   v   - vector of int-s
   n  - dimension of the vector
   r  - number of int-s in one line */
void vec_ffprint_t(FILE *f, double *v, unsigned n, unsigned r) {
  unsigned i;
  for (i=0; i<(n<r ? n : r); i++)
    fprintf(f,"%15d",i+1);
  fprintf(f,"\n");
  for (i=0; i<n; i++) {
    if (((i+1)%r)==1)
      fprintf(f,"%5d",i+1);
    fprintf(f,"%15e",v[i]);
    if (!((i+1)%r))
      fprintf(f,"\n");
    }
  if (n%r)
    fprintf(f,"\n");
  }

/* -------------------------------------------------------------------------- */

/* formatted printing of string vector

   file  - pointer to open file
   vec   - vector of strings
   ntot  - dimension of the vector
   nrow  - number of strings in one line
   width - space for one string */
void vec_sfprint(FILE *file, char **vec, unsigned ntot,
  unsigned nrow, unsigned width) {
  unsigned i;
  for (i=0; i<ntot; i++) {
    fprintf(file,"%.*s",width,vec[i]);
    if ((i+1)%nrow==0)
      fprintf(file,"\n");
    }
  if (ntot%nrow)
    fprintf(file,"\n");
  }

/* -------------------------------------------------------------------------- */
