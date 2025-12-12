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
#include <stdio.h>
#include "cmn/matrix.h"

/* -------------------------------------------------------------------------- */

/* print integer matrix

   mat  - integer matrix
   nrow - number of rows of the matrix
   ncol - number of columns of the matrix */
void mat_iprint(int **mat, unsigned nrow, unsigned ncol) {
  unsigned i,j;
  for (i=0; i<nrow; i++) {
    for (j=0; j<ncol; j++)
      printf("%10d",mat[i][j]);
    printf("\n");
    }
  }

/* print unsigned integer matrix

   mat  - unsigned integer matrix
   nrow - number of rows of the matrix
   ncol - number of columns of the matrix */
void mat_uprint(unsigned **mat, unsigned nrow, unsigned ncol) {
  unsigned i,j;
  for (i=0; i<nrow; i++) {
    for (j=0; j<ncol; j++)
      printf("%10u",mat[i][j]);
    printf("\n");
    }
  }

/* print real matrix

   mat  - integer matrix
   nrow - number of rows of the matrix
   ncol - number of columns of the matrix */
void mat_fprint(double **mat, unsigned nrow, unsigned ncol) {
  unsigned i,j;
  for (i=0; i<nrow; i++) {
    for (j=0; j<ncol; j++)
      printf("%15e",mat[i][j]);
    printf("\n");
    }
  }

/* print real matrix to specific number of columns

   mat  - symmetric double precision real matrix
   nrow - number of rows of the matrix
   ncol - number of columns of the matrix
   nval - number of printed values in line */
void mat_fprint_c(double **mat, unsigned nrow, unsigned ncol, unsigned nval) {
  mat_ffprint_c(stdout,mat,nrow,ncol,nval);
  }

/* print real matrix to specific number of columns to file

   f    - open file stream
   mat  - symmetric double precision real matrix
   nrow - number of rows of the matrix
   ncol - number of columns of the matrix
   nval - number of printed values in line */
void mat_ffprint_c(FILE *f, double **mat, unsigned nrow, 
  unsigned ncol, unsigned nval) {
  mat_ffprint_c_scale(f,mat,nrow,ncol,nval,1.0);
  }

/* print scaled real matrix to specific number of columns to file

   f    - open file stream
   mat  - symmetric double precision real matrix
   nrow - number of rows of the matrix
   ncol - number of columns of the matrix
   nval - number of printed values in line
   scl  - the scaling factor  */
void mat_ffprint_c_scale(FILE *f, double **mat, unsigned nrow, 
  unsigned ncol, unsigned nval, double scl) {
  unsigned i,j,m=0,num=0;
  while (num<ncol) {
     if (num+nval<ncol)
        m = num+nval;
     else
        m = ncol;
     for (i=num; i<m; i++)
        fprintf(f,"%15i",i+1);
     fprintf(f,"\n");
     for (i=0; i<nrow; i++) {
        fprintf(f,"%5i",i+1);
        for (j=num; j<m; j++)
           fprintf(f,"%15e",scl*mat[i][j]);
        fprintf(f,"\n");
        }
     num += nval;
     }
  }

/* print transposed real matrix to specific number of columns

   mat  - symmetric double precision real matrix
   nrow - number of rows of the matrix
   ncol - number of columns of the matrix
   nval - number of printed values in line */
void mat_fprint_tc(double **mat, unsigned nrow, unsigned ncol, unsigned nval) {
  unsigned i,j,m=0,num=0;
  while (num<nrow) {
     if (num+nval<nrow)
        m = num+nval;
     else
        m = nrow;
     for (i=num; i<m; i++)
        printf("%15i",i+1);
     printf("\n");
     for (i=0; i<ncol; i++) {
        printf("%5i",i+1);
        for (j=num; j<m; j++)
           printf("%15e",mat[j][i]);
        printf("\n");
        }
     num += nval;
     }
  }

/* print double-precision real matrix with specific number format
    
   m     - the matrix
   n1,n2 - dimensions of the matrix
   f     - the number format */
void mat_fprint_f(double **m, unsigned n1, unsigned n2, char *f) {
  unsigned i,j;
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++)
      printf(f,m[i][j]);
    printf("\n");
    }
  }

/* print double precision real sub-matrix

   mat     - double precision real matrix
   nr1,nr2 - first and last row for printing
   nc1,nc2 - first and last column for printing */
void mat_fprint_sub(double **mat, unsigned nr1, unsigned nc1,
  unsigned nr2, unsigned nc2) {
  unsigned i,j;
  for (i=nr1; i<nr2; i++) {
    for (j=nc1; j<nc2; j++)
      printf("%15e",mat[i][j]);
    printf("\n");
    }
  }

/* print square real matrix to specific number of columns

   mat - square double precision real matrix
   n   - dimension of the matrix
   k   - number of printed values in line */
void mat_fprint_sc(double **mat, unsigned n, unsigned k) {
  unsigned i,j,m=0,num=0;
  while (num<n) {
     if (num+k<n)
        m = num+k;
     else
        m = n;
     for (i=num+1; i<=m; i++)
        printf("%15i",i);
     printf("\n");
     for (i=0; i<n; i++) {
        printf("%5i",i+1);
        for (j=num; j<m; j++)
           printf("%15e",mat[i][j]);
        printf("\n");
        }
     num += k;
     }
  }

/* print real square matrix to specific number of columns into a file

   file - pointer to open file
   mat  - square double precision real matrix
   n    - dimension of the matrix
   k    - number of printed values in line */
void mat_ffprint_sc(FILE *file, double **mat, unsigned n, unsigned k) {
  unsigned i,j,m=0,num=0;
  while (num<n) {
     if (num+k<n)
        m = num+k;
     else
        m = n;
     for (i=num+1; i<=m; i++)
        fprintf(file,"%15i",i);
     fprintf(file,"\n");
     for (i=0; i<n; i++) {
        fprintf(file,"%5i",i+1);
        for (j=num; j<m; j++)
           fprintf(file,"%15e",mat[i][j]);
        fprintf(file,"\n");
        }
     num += k;
     }
  }

/* print real symmetric matrix to specific number of columns

   mat - symmetric double precision real matrix
   n   - dimension of the matrix
   k   - number of printed values in line */
void mat_fprint_Sc(double **mat, unsigned n, unsigned k) {
  unsigned i,j,m=0,num=0;
  while (num<n) {
     if (num+k<n)
        m = num+k;
     else
        m = n;
     for (i=num; i<m; i++)
        printf("%15i",i+1);
     printf("\n");
     for (i=num; i<n; i++) {
        printf("%5i",i+1);
        if (i<num+k)
           m = i+1;
        else
           m = num+k;
        for (j=num; j<m; j++)
           printf("%15e",mat[i][j]);
        printf("\n");
        }
     num += k;
     }
  }

/* print real symmetric matrix to specific number of columns into a file

   file - pointer to open file
   mat  - symmetric double precision real matrix
   n    - dimension of the matrix
   k    - number of printed values in line */
void mat_ffprint_Sc(FILE *file, double **mat, unsigned n, unsigned k) {
  unsigned i,j,m=0,num=0;
  while (num<n) {
     if (num+k<n)
        m = num+k;
     else
        m = n;
     for (i=num; i<m; i++)
        fprintf(file,"%15i",i+1);
     fprintf(file,"\n");
     for (i=num; i<n; i++) {
        fprintf(file,"%5i",i+1);
        if (i<num+k)
           m = i+1;
        else
           m = num+k;
        for (j=num; j<m; j++)
           fprintf(file,"%15e",mat[i][j]);
        fprintf(file,"\n");
        }
     num += k;
     }
  }

/* print scaled real symmetric matrix to specific number of columns into a file

   file - pointer to open file
   mat  - symmetric double precision real matrix
   n    - dimension of the matrix
   k    - number of printed values in line
   scale - scaling factor */
void mat_ffprint_Sc_scale(FILE *file, double **mat, unsigned n, unsigned k,
  double scale) {
  unsigned i,j,m=0,num=0;
  while (num<n) {
     if (num+k<n)
        m = num+k;
     else
        m = n;
     for (i=num; i<m; i++)
        fprintf(file,"%15i",i+1);
     fprintf(file,"\n");
     for (i=num; i<n; i++) {
        fprintf(file,"%5i",i+1);
        if (i<num+k)
           m = i+1;
        else
           m = num+k;
        for (j=num; j<m; j++)
           fprintf(file,"%15e",mat[i][j]*scale);
        fprintf(file,"\n");
        }
     num += k;
     }
  }

/* -------------------------------------------------------------------------- */

/* print complex matrix

   m  - integer matrix
   nr - number of rows of the matrix
   nc - number of columns of the matrix */
void mat_zprint(complex double **m, unsigned nr, unsigned nc) {
  unsigned i,j;
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) {
      if (cimag(m[i][j])<0.0)
        printf("%15e - i%13e",creal(m[i][j]),fabs(cimag(m[i][j])));
      else
        printf("%15e + i%13e",creal(m[i][j]),fabs(cimag(m[i][j])));
      }
    printf("\n");
    }
  }

/* print complex matrix to specific number of columns

   mat  - symmetric double precision real matrix
   nrow - number of rows of the matrix
   ncol - number of columns of the matrix
   nval - number of printed values in line */
void mat_zprint_c(complex double **mat, unsigned nrow, unsigned ncol,
  unsigned nval) {
  mat_fzprint_c(stdout,mat,nrow,ncol,nval);
  }

/* print complex matrix with specific number format
    
   m     - the matrix
   n1,n2 - dimensions of the matrix
   f     - the number format */
void mat_zprint_f(complex double **m, unsigned n1, unsigned n2, char *f) {
  unsigned i,j;
  char fn[80],fp[80];
  sprintf(fp,"%s + i%s",f,f);
  sprintf(fn,"%s - i%s",f,f);
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++)
      if (cimag(m[i][j])<0.0)
        printf(fn,creal(m[i][j]),fabs(cimag(m[i][j])));
      else
        printf(fp,creal(m[i][j]),fabs(cimag(m[i][j])));
    printf("\n");
    }
  }

/* print complex matrix to specific number of columns to file

   f    - open file stream
   mat  - symmetric double precision real matrix
   nrow - number of rows of the matrix
   ncol - number of columns of the matrix
   nval - number of printed values in line */
void mat_fzprint_c(FILE *f, complex double **mat, unsigned nrow, 
  unsigned ncol, unsigned nval) {
  mat_fzprint_c_scale(f,mat,nrow,ncol,nval,1.0);
  }

/* print scaled complex matrix to specific number of columns to file

   f    - open file stream
   mat  - symmetric double precision real matrix
   nrow - number of rows of the matrix
   ncol - number of columns of the matrix
   nval - number of printed values in line
   scl  - the scaling factor  */
void mat_fzprint_c_scale(FILE *f, complex double **mat, unsigned nrow, 
  unsigned ncol, unsigned nval, complex double scl) {
  unsigned i,j,m=0,num=0;
  complex double t;
  while (num<ncol) {
    if (num+nval<ncol)
      m = num+nval;
    else
      m = ncol;
    for (i=num; i<m; i++)
      fprintf(f,i>0 ? "%25i" : "%18i",i+1);
    fprintf(f,"\n");
    for (i=0; i<nrow; i++) {
       fprintf(f,"%5i",i+1);
       for (j=num; j<m; j++) {
         t = scl*mat[i][j];
         fprintf(f,"  %10.3e %+10.3e i",creal(t),cimag(t));
         }
       fprintf(f,"\n");
       }
    num += nval;
    }
  }

/* -------------------------------------------------------------------------- */
