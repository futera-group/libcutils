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

#include <stdlib.h>
#include "cmn/message.h"
#include "cmn/types.h"

/* -------------------------------------------------------------------------- */

/* allocate rank-1 data array in config item struct

   d    - data pointer
   type - type of data array 
   n    - first data-array dimension */
void cfg_alloc_rank1(void *d, short type, unsigned n) {
  unsigned i;
  int *ti = NULL;
  unsigned *tu = NULL;
  short *th = NULL;
  unsigned short *tg = NULL;
  long *tl = NULL;
  unsigned long *tm = NULL;
  double *td = NULL;
  long double *tx = NULL;
  char **ts = NULL;
  switch (type) {
    case TYPE_INT: /* integer */
      ti = (int*)malloc(n*sizeof(int));
      if (!ti) msg_error("cannot allocate rank-1 int array",1);
      for (i=0; i<n; i++) ti[i] = 0;
      (*((int**)d)) = ti;
      break;
    case TYPE_UINT: /* unsigned integer */
      tu = (unsigned*)malloc(n*sizeof(unsigned));
      if (!tu) msg_error("cannot allocate rank-1 u-int array",1);
      for (i=0; i<n; i++) tu[i] = 0;
      (*((unsigned**)d)) = tu;
      break;
    case TYPE_SINT: /* short integer */
      th = (short*)malloc(n*sizeof(short));
      if (!th) msg_error("cannot allocate rank-1 s-int array",1);
      for (i=0; i<n; i++) th[i] = 0;
      (*((short**)d)) = th;
      break;
    case TYPE_USINT: /* unsigned short integer */
      tg = (unsigned short*)malloc(n*sizeof(unsigned short));
      if (!tg) msg_error("cannot allocate rank-1 u-s-int array",1);
      for (i=0; i<n; i++) tg[i] = 0;
      (*((unsigned short**)d)) = tg;
      break;
    case TYPE_LINT: /* long integer */
      tl = (long*)malloc(n*sizeof(long));
      if (!tl) msg_error("cannot allocate rank-1 l-int array",1);
      for (i=0; i<n; i++) tl[i] = 0;
      (*((long**)d)) = tl;
      break;
    case TYPE_ULINT: /* unsigned long integer */
      tm = (unsigned long*)malloc(n*sizeof(unsigned long));
      if (!tm) msg_error("cannot allocate rank-1 u-l-int array",1);
      for (i=0; i<n; i++) tm[i] = 0;
      (*((unsigned long**)d)) = tm;
      break;
    case TYPE_DOUBLE: /* double-precision real */
      td = (double*)malloc(n*sizeof(double));
      if (!td) msg_error("cannot allocate rank-1 double array",1);
      for (i=0; i<n; i++) td[i] = 0.0;
      (*((double**)d)) = td;
      break;
    case TYPE_LDOUBLE: /* long double-precision real */
      tx = (long double*)malloc(n*sizeof(long double));
      if (!tx) msg_error("cannot allocate rank-1 long double array",1);
      for (i=0; i<n; i++) tx[i] = 0.0;
      (*((long double**)d)) = tx;
      break;
    case TYPE_STRING: /* string */
      ts = (char**)malloc(n*sizeof(char*));
      if (!ts) msg_error("cannot allocate rank-1 string array",1);
      for (i=0; i<n; i++) ts[i] = NULL;
      (*((char***)d)) = ts;
      break;
    default:
      msg_error_f("unsupported config-file rank-1 data (%d)",1,type);
    }
  }

/* allocate rank-2 data array in config item struct

   d    - data pointer
   type - type of data array 
   n    - first data-array dimension */
void cfg_alloc_rank2(void *d, short type, unsigned n) {
  unsigned i;
  int **ti = NULL;
  unsigned **tu = NULL;
  short **th = NULL;
  unsigned short **tg = NULL;
  long **tl = NULL;
  unsigned long **tm = NULL;
  double **td = NULL;
  long double **tx = NULL;
  char ***ts = NULL;
  switch (type) {
    case TYPE_INT: /* integer */
      ti = (int**)malloc(n*sizeof(int*));
      if (!ti) msg_error("cannot allocate rank-2 int array",1);
      for (i=0; i<n; i++) ti[i] = NULL;
      (*((int***)d)) = ti;
      break;
    case TYPE_UINT: /* unsigned integer */
      tu = (unsigned**)malloc(n*sizeof(unsigned*));
      if (!tu) msg_error("cannot allocate rank-2 u-int array",1);
      for (i=0; i<n; i++) tu[i] = NULL;
      (*((unsigned***)d)) = tu;
      break;
    case TYPE_SINT: /* short integer */
      th = (short**)malloc(n*sizeof(short*));
      if (!th) msg_error("cannot allocate rank-2 s-int array",1);
      for (i=0; i<n; i++) th[i] = NULL;
      (*((short***)d)) = th;
      break;
    case TYPE_USINT: /* unsigned short integer */
      tg = (unsigned short**)malloc(n*sizeof(unsigned short*));
      if (!tg) msg_error("cannot allocate rank-2 u-s-int array",1);
      for (i=0; i<n; i++) tg[i] = NULL;
      (*((unsigned short***)d)) = tg;
      break;
    case TYPE_LINT: /* long integer */
      tl = (long**)malloc(n*sizeof(long*));
      if (!tl) msg_error("cannot allocate rank-2 l-int array",1);
      for (i=0; i<n; i++) tl[i] = NULL;
      (*((long***)d)) = tl;
      break;
    case TYPE_ULINT: /* unsigned long integer */
      tm = (unsigned long**)malloc(n*sizeof(unsigned long*));
      if (!tm) msg_error("cannot allocate rank-2 u-l-int array",1);
      for (i=0; i<n; i++) tm[i] = NULL;
      (*((unsigned long***)d)) = tm;
      break;
    case TYPE_DOUBLE: /* double-precision real */
      td = (double**)malloc(n*sizeof(double*));
      if (!td) msg_error("cannot allocate rank-2 double array",1);
      for (i=0; i<n; i++) td[i] = NULL;
      (*((double***)d)) = td;
      break;
    case TYPE_LDOUBLE: /* long double-precision real */
      tx = (long double**)malloc(n*sizeof(long double*));
      if (!tx) msg_error("cannot allocate rank-2 double array",1);
      for (i=0; i<n; i++) tx[i] = NULL;
      (*((long double***)d)) = tx;
      break;
    case TYPE_STRING: /* string */
      ts = (char***)malloc(n*sizeof(char**));
      if (!ts) msg_error("cannot allocate rank-2 string array",1);
      for (i=0; i<n; i++) ts[i] = NULL;
      (*((char****)d)) = ts;
      break;
    default:
      msg_error_f("unsupported config-file rank-2 data (%d)",1,type);
    }
  }

/* allocate rank-3 data array in config item struct

   d    - data pointer
   type - type of data array 
   n    - first data-array dimension */
void cfg_alloc_rank3(void *d, short type, unsigned n) {
  unsigned i;
  int ***ti = NULL;
  unsigned ***tu = NULL;
  short ***th = NULL;
  unsigned short ***tg = NULL;
  long ***tl = NULL;
  unsigned long ***tm = NULL;
  double ***td = NULL;
  long double ***tx = NULL;
  char ****ts = NULL;
  switch (type) {
    case TYPE_INT: /* integer */
      ti = (int***)malloc(n*sizeof(int**));
      if (!ti) msg_error("cannot allocate rank-3 int array",1);
      for (i=0; i<n; i++) ti[i] = NULL;
      (*((int****)d)) = ti;
      break;
    case TYPE_UINT: /* unsigned integer */
      tu = (unsigned***)malloc(n*sizeof(unsigned**));
      if (!tu) msg_error("cannot allocate rank-3 u-int array",1);
      for (i=0; i<n; i++) tu[i] = NULL;
      (*((unsigned****)d)) = tu;
      break;
    case TYPE_SINT: /* short integer */
      th = (short***)malloc(n*sizeof(short**));
      if (!th) msg_error("cannot allocate rank-3 s-int array",1);
      for (i=0; i<n; i++) th[i] = NULL;
      (*((short****)d)) = th;
      break;
    case TYPE_USINT: /* unsigned short integer */
      tg = (unsigned short***)malloc(n*sizeof(unsigned short**));
      if (!tg) msg_error("cannot allocate rank-3 u-s-int array",1);
      for (i=0; i<n; i++) tg[i] = NULL;
      (*((unsigned short****)d)) = tg;
      break;
    case TYPE_LINT: /* long integer */
      tl = (long***)malloc(n*sizeof(long**));
      if (!tl) msg_error("cannot allocate rank-3 l-int array",1);
      for (i=0; i<n; i++) tl[i] = NULL;
      (*((long****)d)) = tl;
      break;
    case TYPE_ULINT: /* unsigned long integer */
      tm = (unsigned long***)malloc(n*sizeof(unsigned long**));
      if (!tm) msg_error("cannot allocate rank-3 u-l-int array",1);
      for (i=0; i<n; i++) tm[i] = NULL;
      (*((unsigned long****)d)) = tm;
      break;
    case TYPE_DOUBLE: /* double-precision real */
      td = (double***)malloc(n*sizeof(double**));
      if (!td) msg_error("cannot allocate rank-3 double array",1);
      for (i=0; i<n; i++) td[i] = NULL;
      (*((double****)d)) = td;
      break;
    case TYPE_LDOUBLE: /* long double-precision real */
      tx = (long double***)malloc(n*sizeof(long double**));
      if (!tx) msg_error("cannot allocate rank-3 long double array",1);
      for (i=0; i<n; i++) tx[i] = NULL;
      (*((long double****)d)) = tx;
      break;
    case TYPE_STRING: /* string */
      ts = (char****)malloc(n*sizeof(char***));
      if (!ts) msg_error("cannot allocate rank-3 string array",1);
      for (i=0; i<n; i++) ts[i] = NULL;
      (*((char*****)d)) = ts;
      break;
    default:
      msg_error_f("unsupported config-file rank-3 data (%d)",1,type);
    }
  }

/* allocate rank-4 data array in config item struct

   d    - data pointer
   type - type of data array 
   n    - first data-array dimension */
void cfg_alloc_rank4(void *d, short type, unsigned n) {
  unsigned i;
  int ****ti = NULL;
  unsigned ****tu = NULL;
  short ****th = NULL;
  unsigned short ****tg = NULL;
  long ****tl = NULL;
  unsigned long ****tm = NULL;
  double ****td = NULL;
  long double ****tx = NULL;
  char *****ts = NULL;
  switch (type) {
    case TYPE_INT: /* integer */
      ti = (int****)malloc(n*sizeof(int***));
      if (!ti) msg_error("cannot allocate rank-4 int array",1);
      for (i=0; i<n; i++) ti[i] = NULL;
      (*((int*****)d)) = ti;
      break;
    case TYPE_UINT: /* unsigned integer */
      tu = (unsigned****)malloc(n*sizeof(unsigned***));
      if (!tu) msg_error("cannot allocate rank-4 u-int array",1);
      for (i=0; i<n; i++) tu[i] = NULL;
      (*((unsigned*****)d)) = tu;
      break;
    case TYPE_SINT: /* short integer */
      th = (short****)malloc(n*sizeof(short***));
      if (!th) msg_error("cannot allocate rank-4 s-int array",1);
      for (i=0; i<n; i++) th[i] = NULL;
      (*((short*****)d)) = th;
      break;
    case TYPE_USINT: /* unsigned short integer */
      tg = (unsigned short****)malloc(n*sizeof(unsigned short***));
      if (!tg) msg_error("cannot allocate rank-4 u-s-int array",1);
      for (i=0; i<n; i++) tg[i] = NULL;
      (*((unsigned short*****)d)) = tg;
      break;
    case TYPE_LINT: /* long integer */
      tl = (long****)malloc(n*sizeof(long***));
      if (!tl) msg_error("cannot allocate rank-4 l-int array",1);
      for (i=0; i<n; i++) tl[i] = NULL;
      (*((long*****)d)) = tl;
      break;
    case TYPE_ULINT: /* unsigned long integer */
      tm = (unsigned long****)malloc(n*sizeof(unsigned long***));
      if (!tm) msg_error("cannot allocate rank-4 u-l-int array",1);
      for (i=0; i<n; i++) tm[i] = NULL;
      (*((unsigned long*****)d)) = tm;
      break;
    case TYPE_DOUBLE: /* double-precision real */
      td = (double****)malloc(n*sizeof(double***));
      if (!td) msg_error("cannot allocate rank-4 double array",1);
      for (i=0; i<n; i++) td[i] = NULL;
      (*((double*****)d)) = td;
      break;
    case TYPE_LDOUBLE: /* long double-precision real */
      tx = (long double****)malloc(n*sizeof(long double***));
      if (!tx) msg_error("cannot allocate rank-4 long double array",1);
      for (i=0; i<n; i++) tx[i] = NULL;
      (*((long double*****)d)) = tx;
      break;
    case TYPE_STRING: /* string */
      ts = (char*****)malloc(n*sizeof(char****));
      if (!ts) msg_error("cannot allocate rank-4 string array",1);
      for (i=0; i<n; i++) ts[i] = NULL;
      (*((char******)d)) = ts;
      break;
    default:
      msg_error_f("unsupported config-file rank-4 data (%d)",1,type);
    }
  }

/* allocate data array in config item struct

   d    - data pointer
   type - type of data array 
   rank - dimension of the array
   n    - first data-array dimension */
void cfg_alloc_data(void *d, short type, int rank, unsigned n) {
  if (d && rank>0) {
    switch (rank) {
      case 1: cfg_alloc_rank1(d,type,n); break;
      case 2: cfg_alloc_rank2(d,type,n); break;
      case 3: cfg_alloc_rank3(d,type,n); break;
      case 4: cfg_alloc_rank4(d,type,n); break;
      default:
        msg_error_f("attempt to allocate config-file data of"
         " unsupported rank (%d)",1,rank);
      }
    }
  }

/* -------------------------------------------------------------------------- */
