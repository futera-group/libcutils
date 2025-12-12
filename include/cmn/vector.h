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

#ifndef ZF_LIB_CMN_VECTOR_H
#define ZF_LIB_CMN_VECTOR_H

#include <complex.h>
#include <stdio.h>

/* -------------------------------------------------------------------------- */

/* symbolic constants */
#define VEC_FRM_f 1
#define VEC_FRM_F 2
#define VEC_FRM_e 3
#define VEC_FRM_E 4
#define VEC_FRM_g 5
#define VEC_FRM_G 6

/* -------------------------------------------------------------------------- */

/* functions for allocation and freeing memory */

/* allocate memory for vector of int-s */
int *vec_ialloc(long unsigned);
/* allocate memory for vector of unsigned int-s */
unsigned *vec_ualloc(long unsigned);
/* allocate memory for vector of short int-s */
short *vec_sialloc(long unsigned);
/* allocate memory for vector of short unsigned int-s */
short unsigned *vec_sualloc(long unsigned);
/* allocate memory for vector of long integers */
long *vec_lialloc(long unsigned);
/* allocate memory for vector of unsigned longs */
unsigned long *vec_lualloc(long unsigned);
/* allocate memory for vector of double precision reals */
double *vec_falloc(long unsigned);
/* allocate memory for vector of longd double precision reals */
long double *vec_lfalloc(long unsigned);
/* allocate memory for vector of double precision complex numbers */
double complex *vec_zalloc(long unsigned);
/* allocate memory for array of strings */
char **vec_salloc(long unsigned, long unsigned);
/* allocate memory for array of user defined structs */
void *vec_talloc(long unsigned, long unsigned);

/* deallocate memory of vector of int-s */
void* vec_ifree(int*);
/* deallocate memory of vector of unsigned int-s */
void* vec_ufree(unsigned*);
/* deallocate memory of vector of short int-s */
void* vec_sifree(short*);
/* deallocate memory of vector of short unsigned int-s */
void* vec_sufree(short unsigned*);
/* deallocate memory of vector of long integers */
void* vec_lifree(long*);
/* deallocate memory of vector of unsigned longs */
void* vec_lufree(unsigned long*);
/* deallocate memory of vector of double precision reals */
void* vec_ffree(double*);
/* deallocate memory of vector of long double precision reals */
void* vec_lffree(long double*);
/* deallocate memory of vector of double precision complex numbers */
void* vec_zfree(double complex*);
/* deallocate memory of array of strings */
void* vec_sfree(char**, long unsigned);
/* deallocate memory of vector of user defined structs */
void* vec_tfree(void*);


/* functions for copying */

/* copy values from one vector to another (integer version) */
void vec_icopy(int*, int*, unsigned);
/* create new copy of vector (integer version) */
int* vec_icopy_new(int*, unsigned);

/* copy values from one vector to another (unsigned integer version) */
void vec_ucopy(unsigned*, unsigned*, unsigned);
/* create new copy of vector (unsigned integer version) */
unsigned* vec_ucopy_new(unsigned*, unsigned);

/* copy values from one vector to another (short integer version) */
void vec_sicopy(short*, short*, unsigned);
/* create new copy of vector (short integer version) */
short* vec_sicopy_new(short*, unsigned);

/* copy values from one vector to another (unsigned long integer version) */
void vec_lucopy(unsigned long*, unsigned long*, unsigned);
/* create new copy of vector (unsigned long integer version) */
unsigned long* vec_lucopy_new(unsigned long*, unsigned);

/* copy values from one vector to another (double precision real version) */
void vec_fcopy(double*, double*, unsigned);
/* copy values from one vector to another (long double precision rea) */
void vec_lfcopy(long double*, long double*, unsigned);
/* create new copy of vector (double precision version) */
double* vec_fcopy_new(double*, unsigned);
/* copy values from one vector to another and scale them (double) */
void vec_fcopy_scaled(double*, double*, double, unsigned);

/* copy values from one vector to another (double precision complex version) */
void vec_zcopy(double complex*, double complex*, unsigned);
/* create new vector and fill it with values from another one (complex) */
double complex* vec_zcopy_new(double complex*, unsigned);

/* copy values from one vector to another (string version) */
void vec_scopy(char**, char**, unsigned);
/* create new copy of vector (string version) */
char** vec_scopy_new(char**, unsigned, unsigned);

/* copy values from one vector to another (user defined type version) */
void vec_tcopy(void*, void*, void(f)(void*,void*), unsigned, unsigned);
/* create new copy of vector (user defined type version) */
void* vec_tcopy_new(void*, void(f)(void*,void*), unsigned, unsigned);

/* functions for vector operations */

/* set the whole integer vector on one value */
void vec_iset(int*, int, long unsigned);
/* set the whole unsigned integer vector on one value */
void vec_uset(unsigned*, unsigned, long unsigned);
/* set the whole short integer vector on one value */
void vec_siset(short*, short, long unsigned);
/* set the whole long integer vector on one value */
void vec_liset(long*, long, long unsigned);
/* set the whole unsigned long integer vector on one value */
void vec_luset(unsigned long*, unsigned long, long unsigned);
/* set the whole floating point vector on one value */
void vec_fset(double*, double, long unsigned);
/* set the whole floating point vector on one value */
void vec_lfset(long double*, long double, long unsigned);
/* set the whole complex real vector on one value */
void vec_zset(complex double*, complex double, long unsigned);
/* set the whole string vector on one value */
void vec_sset(char**, char*, long unsigned);

/* search given value in vector (integer version) */
int vec_ifind(int*, int, unsigned);
/* search given value in vector (unsigned integer version) */
int vec_ufind(unsigned*, unsigned, unsigned);
/* search given value in vector (short integer version) */
int vec_sifind(short*, short, unsigned);
/* search given value in vector (unsigned long integer version) */
int vec_luin(unsigned long*, unsigned long, unsigned);
/* search given value in vector (real version) */
int vec_ffind(double*, double, unsigned);
/* search given value in vector (string version) */
int vec_sfind(char**, char*, unsigned);
/* search given value in vector (user defined types) */
int vec_tfind(void*, void*, int(f)(void*,void*), unsigned, unsigned);

/* shift vector elements right (unsigned int) */
void vec_urshift(unsigned*, unsigned);
/* insert value and shift vector elements right (unsigned int) */
void vec_urshift_id(unsigned*, unsigned, unsigned, unsigned);
/* shift vector elements left (double precision real) */
void vec_flshift(double*, unsigned);
/* shift vector elements right (double precision real) */
void vec_frshift(double*, unsigned);
/* insert value and shift vector elements right (double) */
void vec_frshift_id(double*, unsigned, unsigned, double);

/* find the smallest value in given vector (unsinged int version) */
unsigned vec_umin(unsigned*, unsigned);
/* find the largest value in given vector (unsinged int version) */
unsigned vec_umax(unsigned*, unsigned);
/* find minimum value in given vector (double precision real) */
unsigned vec_fmin(double*, unsigned);
/* find minimum abs value in given vector (double precision real) */
unsigned vec_fmin_abs(double*, unsigned);
/* find minimum and maximum from values in the vector (double) */
void vec_fminmax(double*, long unsigned, double*, double*);
/* find maximum value in given vector (double precision real) */
unsigned vec_fmax(double*, unsigned);
/* find maximum abs value in given vector (double precision real) */
unsigned vec_fmax_abs(double*, unsigned);
/* find longest string  in given vector (string version) */
unsigned vec_smax(char**, unsigned);

/* norm of the floating point vector */
double vec_fnorm(double*, unsigned);
/* normalize vector to unit (double precision real version) */
void vec_funit(double*, unsigned);
/* create new normalized vector form old one (double precision real version) */
void vec_funit_copy(double*, double*, unsigned);

/* scalar product of two floating point vectors */
double vec_fprod_scalar(double*, double*, unsigned);
/* vector product of two 3D floating point vectors */
void vec_fprod_vector(double*, double*, double*);
/* mixed product of three 3D floating point vectors */
double vec_fprod_mixed(double*, double*, double*);

/* angle between two floating point vectors */
double vec_fangle(double*, double*, unsigned);

/* sum two vectors together (double precision real version) */
void vec_fadd(double*, double*, double*, unsigned);
/* scale and add one vector to another (double precision real version) */
void vec_fadd_scaled(double*, double*, double, unsigned);
/* subtract two vectors (double precision real version) */
void vec_fsub(double*, double*, double*, unsigned);
/* return sum of all elements in the vector (double precision real version) */
double vec_fsum(double*, unsigned);
/* return sum of absolute values of all elements (double real version) */
double vec_fsum_abs(double*, unsigned);
/* scale all vector elements by given factor (double precision real version) */
void vec_fscale(double*, double, unsigned);

/* merge two unsigned integer vector to the first one */
void vec_umerge(unsigned**, unsigned*, unsigned*, unsigned);
/* merge two unsigned integer vector to new one */
unsigned *vec_umerge_new(unsigned*, unsigned, unsigned*, unsigned);

/* function for sorting */

/* sort vector (u-int, bubble sort) */
void vec_usort(unsigned*, unsigned);
/* sort 2 vectors (u1,u2) according to u1 values (bubble sort) */
void vec_uu2sort(unsigned*, unsigned*, unsigned);
/* sort vector (double, bubble sort) */
void vec_fsort(double*, unsigned);
/* sort vector (double, bubble sort) by absolute values of the numbers */
void vec_fsort_abs(double*, unsigned);
/* sort 2 vectors (f,u) according to f values (bubble sort) */
void vec_fu2sort(double*, unsigned*, unsigned);
/* sort 2 vectors (f1,f2) according to f1 values (bubble sort) */
void vec_ff2sort(double*, double*, unsigned);
/* sort vector (string, bubble sort) */
void vec_ssort(char**, unsigned);
/* sort 3 vectors (d1,d2,d3) according to d1 values (bubble sort) */
void vec_t3sort(void*, void*, void*,
  short (*)(void*, unsigned, unsigned),
  void (*)(void*, void*, void*, unsigned, unsigned), unsigned);
/* sort 4 vectors (d1,d2,d3,d4) according to d1 values (bubble sort) */
void vec_t4sort(void*, void*, void*, void*,
  short (*)(void*, unsigned, unsigned),
  void (*)(void*, void*, void*, void*, unsigned, unsigned), unsigned);

/* function for changing size of vectors */

/* reduce vector according to given mask (integer version) */
int *vec_ireduce_mask(int*, unsigned, unsigned*, unsigned);
/* reduce vector to given range (integer version) */
int *vec_ireduce_range(int*, unsigned, unsigned, unsigned);
/* delete multiple occurence of elements (u-int) */
unsigned *vec_ureduce_unique(unsigned*, unsigned*);
/* reduce vector according to given mask (double precision real version) */
double *vec_freduce_mask(double*, unsigned, unsigned*, unsigned);
/* reduce vector to given range (double precision real version) */
double *vec_freduce_range(double*, unsigned, unsigned, unsigned);
/* reduce vector according to given mask (string version) */
char **vec_sreduce_mask(char**, unsigned, unsigned*, unsigned, unsigned);
/* reduce vector to given range (string version) */
char **vec_sreduce_range(char**, unsigned, unsigned, unsigned, unsigned);

/* change dimension of vector (int) */
int *vec_iresize(int*, unsigned, unsigned);
/* change dimension of vector (u-int) */
unsigned *vec_uresize(unsigned*, unsigned, unsigned);
/* change dimension of vector (s-int) */
short *vec_siresize(short*, unsigned, unsigned);
/* change dimension of vector (l-int) */
long *vec_liresize(long*, unsigned, unsigned);
/* change dimension of vector (lu-int) */
long unsigned *vec_luresize(long unsigned*, unsigned, unsigned);
/* change dimension of vector (double) */
double *vec_fresize(double*, unsigned, unsigned);
/* change dimension of vector (string) */
char **vec_sresize(char**, unsigned, unsigned);
/* change dimension of vector (user defined type) */
void *vec_tresize(void*, unsigned, unsigned, unsigned, void(f)(void*,void*));

/* delete one element in vector (integer) */
int *vec_idelete(int*, unsigned, unsigned);
/* delete one element in vector (u-int) */
unsigned *vec_udelete(unsigned*, unsigned, unsigned);
/* delete one element in vector (double) */
double *vec_fdelete(double*, unsigned, unsigned);
/* delete one element in vector (user defined type) */
void *vec_tdelete(void*, unsigned, unsigned, unsigned, void(f)(void*,void*));

/* functions for printing */

/* raw print of integer vector */
void vec_iprint(FILE*, int*, unsigned);
/* formatted printing of integer vector */
void vec_ifprint(FILE*, int*, unsigned, unsigned, unsigned);
/* raw print of short-integer vector */
void vec_siprint(FILE*, short*, unsigned);
/* formatted printing of short-integer vector */
void vec_sifprint(FILE*, short*, unsigned, unsigned, unsigned);

/* raw print of unsigned integer vector */
void vec_uprint(FILE*, unsigned*, unsigned);
/* formatted printing of integer vector */
void vec_ufprint(FILE*, unsigned*, unsigned, unsigned, unsigned);

/* print double precision real vector */
void vec_fprint(double*, unsigned);
/* print double-precision real vector with specific number format */
void vec_fprint_f(double*, unsigned, char*);
/* print scaled double precision real vector */
void vec_fprint_scale(double*, double, unsigned);
/* print double precision real sub-vector */
void vec_fprint_sub(double*, unsigned, unsigned);
/* formatted printing of floating point vector */
void vec_ffprint(FILE*, double*, unsigned, unsigned, unsigned, unsigned, short);
/* formatted printing of double precision float vector to table format */
void vec_ffprint_t(FILE*, double*, unsigned, unsigned);

/* formatted printing of string vector */
void vec_sfprint(FILE*, char**, unsigned, unsigned, unsigned);

/* -------------------------------------------------------------------------- */

#endif
