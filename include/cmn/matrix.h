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

#ifndef ZF_LIB_CMN_MATRIX_H
#define ZF_LIB_CMN_MATRIX_H

#include <stdio.h>
#include <complex.h>

/* -------------------------------------------------------------------------- */

#define MAT_CHECK_PREC     1.0E-10

#define MAT_DIAG_SORT_ASC  1  /* sort eigenvalues in ascending order */
#define MAT_DIAG_SORT_DSC  2  /* sort eigenvalues in descending order */

/* -------------------------------------------------------------------------- */

/* functions allocating and freeing memory */

/* allocate memory for matrix of int-s */
int **mat_ialloc(unsigned, unsigned);
/* allocate memory for matrix of unsigned int-s */
unsigned **mat_ualloc(unsigned, unsigned);
/* allocate memory for matrix of unsigned long int-s */
long unsigned **mat_lualloc(unsigned, unsigned);
/* allocate memory for matrix of short int-s */
short **mat_sialloc(unsigned, unsigned);
/* allocate memory for matrix of double precision floating points */
double **mat_falloc(unsigned, unsigned);
/* allocate memory for complex matrix */
double complex** mat_zalloc(unsigned, unsigned);

/* free memory of matrix of int-s */
void* mat_ifree(int**, unsigned);
/* free memory of matrix of unsigned int-s */
void* mat_ufree(unsigned**, unsigned);
/* free memory of matrix of unsigned long int-s */
void* mat_lufree(long unsigned**, unsigned);
/* free memory of matrix of short int-s */
void* mat_sifree(short int**, unsigned int);
/* free memory of matrix of doubles */
void* mat_ffree(double**, unsigned);
/* free memory allocated for complex matrix */
double complex** mat_zfree(double complex**, unsigned);

/* matrix resize / delete */

/* change number of rows in unsigned integer matrix */
unsigned **mat_uresize_r(unsigned**, unsigned, unsigned, unsigned);
/* change number of rows in floating-point matrix */
double **mat_fresize_r(double**, unsigned, unsigned, unsigned);

/* delete specific row in floating-point matrix */
unsigned **mat_udelete_r(unsigned**, unsigned, unsigned, unsigned);
/* delete specific row in floating-point matrix */
double **mat_fdelete_r(double**, unsigned, unsigned, unsigned);

/* check matrix properties */

/* check if two double real matrices are the same */
short mat_fis_same(double**, double**, unsigned, unsigned);
/* check if two double real matrices are the same and print differences */
short mat_fis_same_print(double**, double**, unsigned, unsigned);
/* check if double real matrix is symmetric */
short mat_fis_sym(double**, unsigned);
/* check if double real matrix is symmetric and print off-diag differences */
short mat_fis_sym_print(double **, unsigned);
/* check if double-precision real matrix is orthogonal */
short mat_fis_orth(double**, unsigned, unsigned);
/* check if double real matrix is orthogonal */
short mat_fis_orth_Sc(double**, unsigned);
/* check if double real matrix is orthogonal and print all vector products */
short mat_fis_orth_Sc_print(double**, unsigned);
/* check if double real matrix is unit matrix */
short mat_fis_unit(double**, unsigned, double*);
/* check if double real matrix is diagonal */
short mat_fis_diag(double**, unsigned);

/* functions for copying matrices */

/* copy values from one matrix to another (integer version) */
void mat_icopy(int**, int**, unsigned, unsigned);
/* create new copy of the matrix (integer version) */
int **mat_icopy_new(int**, unsigned, unsigned);

/* copy values from one matrix to another (short integer version) */
void mat_sicopy(short int**, short int**, unsigned, unsigned);
/* create new copy of the matrix (short integer version) */
short int **mat_sicopy_new(short int**, unsigned, unsigned);

/* copy values from one matrix to another (unsigned integer version) */
void mat_ucopy(unsigned**, unsigned**, unsigned, unsigned);
/* create new copy of the matrix (unsigned integer version) */
unsigned **mat_ucopy_new(unsigned**, unsigned, unsigned);

/* copy values from one matrix to another (double real version) */
void mat_fcopy(double**, double**, unsigned, unsigned);
/* copy transposed values from one matrix to another (double real version) */
void mat_fcopy_trans(double**, double**, unsigned, unsigned);
/* create new copy of the matrix (double version) */
double **mat_fcopy_new(double**, unsigned, unsigned);
/* create new matrix and fill it with transposed values from another one (double version) */
double **mat_fcopy_trans_new(double **, unsigned, unsigned);
/* copy values from sub-matrix (double real version) */
void mat_fcopy_sub(double**, double**, unsigned, unsigned, unsigned, unsigned);

/* copy values from one matrix to another (complex real version) */
void mat_zcopy(complex double**, complex double**, unsigned, unsigned);

/* matrix operation functions */

/* set all elements of matrix to one value (integer version) */
void mat_iset(int**, int, unsigned, unsigned);
/* set all elements of matrix to one value (unsigned int version) */
void mat_uset(unsigned**, unsigned, unsigned, unsigned);
/* set all elements of matrix to one value (short integer version) */
void mat_siset(short int**, short int, unsigned, unsigned);
/* set all elements of matrix to one value (unsigned long int version) */
void mat_luset(long unsigned**, long unsigned, unsigned, unsigned);
/* set all elements of matrix to one value (double precision real version) */
void mat_fset(double**, double, unsigned, unsigned);
/* set all elements of matrix to one value (complex version) */
void mat_zset(double complex**, double complex, unsigned, unsigned);

/* calculate matrix determinant up to size 3 (double-precision real) */
double mat_fdet(double**, unsigned);
/* calculate matrix determinant up to size 3 (double-precision complex) */
complex double mat_zdet(complex double**, unsigned);
/* calculate matrix determinant by Cramer rule (double-precision real) */
double mat_fdet_cramer(double**, unsigned);
/* calculate matrix determinant by Cramer rule (double-precision complex) */
double complex mat_zdet_cramer(double complex**, unsigned);
/* calculate matrix determinant by LU decomposition (double-precision real) */
double mat_fdet_lu(double**, unsigned, unsigned);
/* calculate matrix determinant by LU decomposition (double-prec complex) */
complex double mat_zdet_lu(complex double**, unsigned, unsigned);

/* calculate inverse matrix up to size 3 (double-precision real) */
void mat_finv(double**, double**, unsigned);
/* calculate inverse matrix up to size 3 (double-precision complex) */
void mat_zinv(complex double**, complex double**, unsigned);
/* calculate inverse matrix by Cramer rule (double-precision real) */
void mat_finv_cramer(double**, double**, unsigned);
/* calculate inverse matrix by Cramer rule (double-precision complex) */
void mat_zinv_cramer(double complex**, double complex**, unsigned);
/* calculate inverse matrix by LU decomposition (double-precision real) */
void mat_finv_lu(double**, double**, unsigned*, unsigned);
/* calculate inverse matrix by LU decomposition (double-precision complex) */
void mat_zinv_lu(complex double**, complex double**, unsigned*, unsigned);

/* set square matrix to unit (integer version) */
void mat_iunit(int**, unsigned);
/* set square matrix to unit (unsigned version) */
void mat_uunit(unsigned**, unsigned);
/* set square matrix to unit (double precision real version) */
void mat_funit(double**, unsigned);

/* transpose integer square matrix */
void mat_itrans(int**, unsigned);
/* transpose unsigned integer square matrix */
void mat_utrans(unsigned**, unsigned);
/* transpose real square matrix */
void mat_ftrans(double**, unsigned);

/* multiplication of two matrices (integer version) */
void mat_imult(int**, int**, int**, unsigned, unsigned, unsigned);
/* multiplication of matrix by vector from left (integer version) */
void mat_imult_lvec(int*, int**, int*, unsigned, unsigned);
/* multiplication of matrix by vector from right (integer version) */
void mat_imult_rvec(int*, int**, int*, unsigned, unsigned);

/* multiplication of two matrices (unsigned integer version) */
void mat_umult(unsigned**, unsigned**, unsigned**, unsigned, unsigned,
  unsigned);
/* multiplication of matrix by vector from left (integer version) */
void mat_umult_lvec(unsigned*, unsigned**, unsigned*, unsigned, unsigned);
/* multiplication of matrix by vector from right (unsigned integer version) */
void mat_umult_rvec(unsigned*, unsigned**, unsigned*, unsigned, unsigned);

/* multiplication of two matrices (double precision real version) */
void mat_fmult(double**, double**, double**, unsigned, unsigned, unsigned);
/* multiplication of matrix by vector from left (double real version) */
void mat_fmult_lvec(double*, double**, double*, unsigned, unsigned);
/* multiplication of matrix by vector from right (double real version) */
void mat_fmult_rvec(double*, double**, double*, unsigned, unsigned);

/* multiply two square matrices in place, product in the first one */
void mat_fmult_sqr(double **, double **, unsigned);
/* multiply two square matrices in place, product in the second one */
void mat_fmult_sql(double **, double **, unsigned);
/* multiply two square 4c-format matrices in place, product in the first one */
void mat_fmult_sqr_4cf(double*, double*);
/* multiply two square 4c-format matrices in place, product in the second one */
void mat_fmult_sql_4cf(double*, double*);
/* multiply square matrix by smaller square matrix block from the right */
void mat_fmult_sqsub(double **, double **, unsigned, unsigned);

/* similarity transformation of square matrix in place */
void mat_fmult_sim(double **, double **, unsigned);

/* multiplication of two matrices (complex version) */
void mat_zmult(double complex **, double complex **, double complex **,
  unsigned, unsigned, unsigned);

/* diagonalization */

/* diagonalize square matrix, calculate eigenvalues and eigenvectors */
void mat_diag(double**, double**, double*, unsigned, short);
/* reduce symmetric matrix to tridiagonal form by Householder method */
void mat_diag_householder(double**, double*, double*, unsigned);
/* diagonalize 3-diagonal matrix by QM algorithm with implicit shifts */
void mat_diag_ql_implicit(double*, double*, double**, unsigned);
/* sort eigenvalues into ascending order by straiht insertion method */
void mat_diag_sort_asc(double*, double**, unsigned);
/* sort eigenvalues into descending order by straiht insertion method */
void mat_diag_sort_dsc(double*, double**, unsigned);

/* decomposition */

/* singular value decomposition of general matrix */
void mat_fdecomp_svd(double**, double**, double**, double*, unsigned, unsigned);
/* single value decomposition - initial reduction to bidiagonal form */
void mat_fdecomp_svd_hauseholder(double**, double*, unsigned, unsigned,
  double*, double*, double*, double*);

/* LU decomposition of a square matrix (double-precision real) */
void mat_fdecomp_lu(double**, unsigned, unsigned*, unsigned*);
/* LU decomposition of a square matrix (double-precision complex) */
void mat_zdecomp_lu(complex double**, unsigned, unsigned*, unsigned*);

/* row / column shifting */

/* shift rows down in matrix (integer version) */
void mat_idshift(int**, unsigned, unsigned);
/* shift rows down in matrix (unsigned integer version) */
void mat_udshift(unsigned**, unsigned, unsigned);
/* shift rows down in matrix (double real version) */
void mat_fdshift(double**, unsigned, unsigned);

/* shift rows up in matrix (integer version) */
void mat_iushift(int**, unsigned, unsigned);
/* shift rows up in matrix (unsigned integer version) */
void mat_uushift(unsigned**, unsigned, unsigned);
/* shift rows up in matrix (double real version) */
void mat_fushift(double**, unsigned, unsigned);

/* shift columns left in matrix (integer version) */
void mat_ilshift(int**, unsigned, unsigned);
/* shift columns left in matrix (unsigned integer version) */
void mat_ulshift(unsigned**, unsigned, unsigned);
/* shift columns left in matrix (double real version) */
void mat_flshift(double**, unsigned, unsigned);

/* shift columns right in matrix (integer version) */
void mat_irshift(int**, unsigned, unsigned);
/* shift columns right in matrix (unsigned integer version) */
void mat_urshift(unsigned**, unsigned, unsigned);
/* shift columns right in matrix (double real version) */
void mat_frshift(double**, unsigned, unsigned);

/* printing functions */

/* print integer matrix */
void mat_iprint(int**, unsigned, unsigned);
/* print unsigned integer matrix */
void mat_uprint(unsigned**, unsigned, unsigned);

/* print real matrix */
void mat_fprint(double**, unsigned, unsigned);
/* print double-precision real matrix with specific number format */
void mat_fprint_f(double**, unsigned, unsigned, char*);
/* print real matrix to specific number of columns */
void mat_fprint_c(double**, unsigned, unsigned, unsigned);
/* print real matrix to specific number of columns to file */
void mat_ffprint_c(FILE*, double**, unsigned, unsigned, unsigned);
/* print scaled real matrix to specific number of columns to file */
void mat_ffprint_c_scale(FILE*, double**, unsigned, unsigned, unsigned, double);
/* print transposed real matrix to specific number of columns */
void mat_fprint_tc(double**, unsigned, unsigned, unsigned);
/* print double precision real sub-matrix */
void mat_fprint_sub(double**, unsigned, unsigned, unsigned, unsigned);
/* print real square matrix to specific number of columns */
void mat_fprint_sc(double**, unsigned, unsigned);
/* print real square matrix to specific number of columns into a file */
void mat_ffprint_sc(FILE*, double**, unsigned, unsigned);
/* print real symmetric matrix to specific number of columns */
void mat_fprint_Sc(double**, unsigned, unsigned);
/* print real symmetric matrix to specific number of columns into a file */
void mat_ffprint_Sc(FILE*, double**, unsigned, unsigned);
/* print real scaled symmetric matrix to specific number of columns in a file */
void mat_ffprint_Sc_scale(FILE*, double**, unsigned, unsigned, double);

/* print complex matrix */
void mat_zprint(complex double**, unsigned, unsigned);
/* print complex matrix with specific number format */
void mat_zprint_f(complex double**m, unsigned, unsigned, char*);
/* print complex matrix to specific number of columns */
void mat_zprint_c(complex double**, unsigned, unsigned, unsigned);
/* print complex matrix to specific number of columns to file */
void mat_fzprint_c(FILE *f, complex double**, unsigned, unsigned, unsigned);
/* print scaled complex matrix to specific number of columns to file */
void mat_fzprint_c_scale(FILE*, complex double**, unsigned, unsigned,
  unsigned, complex double);

/* -------------------------------------------------------------------------- */

#endif
