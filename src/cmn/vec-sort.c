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
#include <string.h>

/* -------------------------------------------------------------------------- */

/* sort vector (u-int, bubble sort)
 
   v - vector of unsigned integers
   n - number of elements in the vector */
void vec_usort(unsigned *v, unsigned n) {
  int i,j,change,val;
  for (i=1; i<n; i++) {
    change = 0;
    for (j=n-1; j>=i; j--)
      if (v[j-1]>v[j]) {
        val = v[j-1];
        v[j-1] = v[j];
        v[j] = val;
        change++;
        }
    if (change==0)
      break;
    }
  }

/* sort 2 vectors (u1,u2) according to u1 values (bubble sort)
 
   v - vector of unsigned integers
   n - number of elements in the vector */
void vec_uu2sort(unsigned *u1, unsigned *u2, unsigned n) {
  unsigned i,j,u1v,u2v;
  short change;
  for (i=1; i<n; i++) {
    change = 0;
    for (j=n-1; j>=i; j--)
      if (u1[j-1]>u1[j]) {
        u1v = u1[j-1];
        u2v = u2[j-1];
        u1[j-1] = u1[j];
        u2[j-1] = u2[j];
        u1[j] = u1v;
        u2[j] = u2v;
        change++;
        }
    if (change==0)
      break;
    }
  }

/* -------------------------------------------------------------------------- */

/* sort vector (double, bubble sort)
 
   v - vector of unsigned integers
   n - number of elements in the vector */
void vec_fsort(double *v, unsigned n) {
  int i,j,change;
  double val;
  for (i=1; i<n; i++) {
    change = 0;
    for (j=n-1; j>=i; j--)
      if (v[j-1]>v[j]) {
        val = v[j-1];
        v[j-1] = v[j];
        v[j] = val;
        change++;
        }
    if (change==0)
      break;
    }
  }

/* sort vector (double, bubble sort) by absolute values of the numbers
 
   v - vector of unsigned integers
   n - number of elements in the vector */
void vec_fsort_abs(double *v, unsigned n) {
  int i,j,change;
  double val;
  for (i=1; i<n; i++) {
    change = 0;
    for (j=n-1; j>=i; j--)
      if (fabs(v[j-1])>fabs(v[j])) {
        val = v[j-1];
        v[j-1] = v[j];
        v[j] = val;
        change++;
        }
    if (change==0)
      break;
    }
  }

/* sort 2 vectors (f,u) according to f values (bubble sort)
 
   v - vector of unsigned integers
   n - number of elements in the vector */
void vec_fu2sort(double *f, unsigned *u, unsigned n) {
  short change;
  unsigned i,j,uv;
  double fv;
  for (i=1; i<n; i++) {
    change = 0;
    for (j=n-1; j>=i; j--)
      if (f[j-1]>f[j]) {
        fv = f[j-1];
        uv = u[j-1];
        f[j-1] = f[j];
        u[j-1] = u[j];
        f[j] = fv;
        u[j] = uv;
        change++;
        }
    if (change==0)
      break;
    }
  }

/* sort 2 vectors (f1,f2) according to f1 values (bubble sort)
 
   v - vector of unsigned integers
   n - number of elements in the vector */
void vec_ff2sort(double *f1, double *f2, unsigned n) {
  short change;
  unsigned i,j;
  double v1,v2;
  for (i=1; i<n; i++) {
    change = 0;
    for (j=n-1; j>=i; j--)
      if (f1[j-1]>f1[j]) {
        v1 = f1[j-1];
        v2 = f2[j-1];
        f1[j-1] = f1[j];
        f2[j-1] = f2[j];
        f1[j] = v1;
        f2[j] = v2;
        change++;
        }
    if (change==0)
      break;
    }
  }

/* -------------------------------------------------------------------------- */

/* sort vector (string, bubble sort)
 
   v - vector of strings
   n - number of elements in the vector */
void vec_ssort(char **v, unsigned n) {
  char *val;
  int i,j,change;
  for (i=1; i<n; i++) {
    change = 0;
    for (j=n-1; j>=i; j--)
      if (strcmp(v[j-1],v[j])>0) {
        val = v[j-1];
        v[j-1] = v[j];
        v[j] = val;
        change++;
        }
    if (change==0)
      break;
    }
  }

/* -------------------------------------------------------------------------- */

/* sort 3 vectors (d1,d2,d3) according to d1 values (bubble sort)
 
   d1,d2,d3 - data vectors
   f1       - function for comparing d1 values 
   f2       - function for switching d1,d2,d3 values
   n        - number of elements in the vectors */
void vec_t3sort(void *d1, void *d2, void *d3,
  short (f1)(void*, unsigned, unsigned),
  void (f2)(void*, void*, void*, unsigned, unsigned), unsigned n) {
  unsigned i,j;
  short change;
  for (i=1; i<n; i++) {
    change = 0;
    for (j=n-1; j>=i; j--) 
      if (f1(d1,j-1,j)) {
        f2(d1,d2,d3,j-1,j);
        change++;
        }
    if (change==0)
      break;
    }
  }

/* sort 4 vectors (d1,d2,d3,d4) according to d1 values (bubble sort)
 
   d1-4 - data vectors
   f1   - function for comparing d1 values 
   f2   - function for switching d1,d2,d3 values
   n    - number of elements in the vectors */
void vec_t4sort(void *d1, void *d2, void *d3, void *d4,
  short (f1)(void*, unsigned, unsigned),
  void (f2)(void*, void*, void*, void*, unsigned, unsigned), unsigned n) {
  unsigned i,j;
  short change;
  for (i=1; i<n; i++) {
    change = 0;
    for (j=n-1; j>=i; j--) 
      if (f1(d1,j-1,j)) {
        f2(d1,d2,d3,d4,j-1,j);
        change++;
        }
    if (change==0)
      break;
    }
  }

/* -------------------------------------------------------------------------- */
