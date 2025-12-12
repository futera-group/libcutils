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

#include "cmn/message.h"
#include "cmn/queue.h"
#include "cmn/string.h"
#include "cmn/types.h"
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* convert string to array of numbers with specified type
 
   s - the string with array
   t - type of the numbers
   d - pointer to data storage place
   n - number of values found in the string */
void str_parse_array(char *s, short t, void *d, unsigned *n) {
  switch (t) {
    case TYPE_INT:    (*((int**)d)) = str_parse_iarray(s,n);             break; 
    case TYPE_UINT:   (*((unsigned**)d)) = str_parse_uarray(s,n);        break;
    case TYPE_SINT:   (*((short**)d)) = str_parse_siarray(s,n);          break;
    case TYPE_USINT:  (*((short unsigned**)d)) = str_parse_suarray(s,n); break;
    case TYPE_LINT:   (*((long**)d)) = str_parse_liarray(s,n);           break;
    case TYPE_ULINT:  (*((long unsigned**)d)) = str_parse_luarray(s,n);  break;
    case TYPE_DOUBLE: (*((double**)d)) = str_parse_farray(s,n);          break;
    case TYPE_LDOUBLE: (*((long double**)d)) = str_parse_lfarray(s,n);   break;
    default:
      msg_error("unsupported type for string array conversion",1);
    }
  }

/* convert string to array of int numbers
 
   s - the string with array
   n - number of values found in the string */
int *str_parse_iarray(char *s, unsigned *n) {
  unsigned i,sv_n;
  char **sv;
  int *v;
  /* split to strings */
  sv = str_split(s,' ',&sv_n);
  /* convert to u-int */
  (*n) = 0;
  v = vec_ialloc(sv_n);
  for (i=0; i<sv_n; i++) {
    if (sscanf(sv[i],"%d",&(v[i]))!=1)
      break;
    (*n)++;
    }
  return(v);
  }

/* convert string to array of u-int numbers
 
   s - the string with array
   n - number of values found in the string */
unsigned *str_parse_uarray(char *s, unsigned *n) {
  unsigned i,sv_n,*v;
  char **sv;
  /* split to strings */
  sv = str_split(s,' ',&sv_n);
  /* convert to u-int */
  (*n) = 0;
  v = vec_ualloc(sv_n);
  for (i=0; i<sv_n; i++) {
    if (sscanf(sv[i],"%u",&(v[i]))!=1)
      break;
    (*n)++;
    }
  return(v);
  }

/* convert string to array of s-int numbers
 
   s - the string with array
   n - number of values found in the string */
short *str_parse_siarray(char *s, unsigned *n) {
  unsigned i,sv_n;
  char **sv;
  short *v;
  /* split to strings */
  sv = str_split(s,' ',&sv_n);
  /* convert to u-int */
  (*n) = 0;
  v = vec_sialloc(sv_n);
  for (i=0; i<sv_n; i++) {
    if (sscanf(sv[i],"%hd",&(v[i]))!=1)
      break;
    (*n)++;
    }
  return(v);
  }

/* convert string to array of su-int numbers
 
   s - the string with array
   n - number of values found in the string */
short unsigned *str_parse_suarray(char *s, unsigned *n) {
  short unsigned *v;
  unsigned i,sv_n;
  char **sv;
  /* split to strings */
  sv = str_split(s,' ',&sv_n);
  /* convert to u-int */
  (*n) = 0;
  v = vec_sualloc(sv_n);
  for (i=0; i<sv_n; i++) {
    if (sscanf(sv[i],"%hu",&(v[i]))!=1)
      break;
    (*n)++;
    }
  return(v);
  }

/* convert string to array of l-int numbers
 
   s - the string with array
   n - number of values found in the string */
long *str_parse_liarray(char *s, unsigned *n) {
  unsigned i,sv_n;
  char **sv;
  long *v;
  /* split to strings */
  sv = str_split(s,' ',&sv_n);
  /* convert to u-int */
  (*n) = 0;
  v = vec_lialloc(sv_n);
  for (i=0; i<sv_n; i++) {
    if (sscanf(sv[i],"%ld",&(v[i]))!=1)
      break;
    (*n)++;
    }
  return(v);
  }

/* convert string to array of lu-int numbers
 
   s - the string with array
   n - number of values found in the string */
long unsigned *str_parse_luarray(char *s, unsigned *n) {
  long unsigned *v;
  unsigned i,sv_n;
  char **sv;
  /* split to strings */
  sv = str_split(s,' ',&sv_n);
  /* convert to u-int */
  (*n) = 0;
  v = vec_lualloc(sv_n);
  for (i=0; i<sv_n; i++) {
    if (sscanf(sv[i],"%lu",&(v[i]))!=1)
      break;
    (*n)++;
    }
  return(v);
  }

/* convert string to array of real numbers
 
   s - the string with array
   n - number of values found in the string */
double *str_parse_farray(char *s, unsigned *n) {
  unsigned i,sv_n;
  double *v;
  char **sv;
  /* split to strings */
  sv = str_split(s,' ',&sv_n);
  /* convert to doubles */
  (*n) = 0;
  v = vec_falloc(sv_n);
  for (i=0; i<sv_n; i++) {
    if (sscanf(sv[i],"%lf",&(v[i]))!=1)
      break;
    (*n)++;
    }
  return(v);
  }

/* convert string to array of real numbers (long double)
 
   s - the string with array
   n - number of values found in the string */
long double *str_parse_lfarray(char *s, unsigned *n) {
  unsigned i,sv_n;
  long double *v;
  char **sv;
  /* split to strings */
  sv = str_split(s,' ',&sv_n);
  /* convert to long doubles */
  (*n) = 0;
  v = vec_lfalloc(sv_n);
  for (i=0; i<sv_n; i++) {
    if (sscanf(sv[i],"%Lf",&(v[i]))!=1)
      break;
    (*n)++;
    }
  return(v);
  }

/* -------------------------------------------------------------------------- */

/* convert one number from string format to integer type
 
   s - string representing the number
   a - start pointer for reading
   b - end pointer for reading */
unsigned str_parse_unum_one(char *s, unsigned a, unsigned b) {
  unsigned i=a,num=0;
  while (i<b && s[i]>='0' && s[i]<='9')
    num = num*10+s[i++]-'0';
  return(num);
  }

/* parse string containing list of numbers separated by commas and dashes
 
   s - the string with array
   t - type of the numbers
   d - pointer to data storage place
   n - number of values found in the string */
void str_parse_num(char *s, short t, void *d, unsigned *n) {
  switch (t) {
    case TYPE_UINT: (*((unsigned**)d)) = str_parse_unum(s,n); break;
    default:
      msg_error("unsupported type for general string array conversion",1);
    }
  }

/* parse part of string between two commas, save numbers into queue
 
   s - string representing unsigned integers
   a - start pointer for reading
   b - end pointer for reading
   q - queue for saving numbers read from string */
void str_parse_unum_commas(char *s, unsigned a, unsigned b, struct queue *q) {
  unsigned i,j,n1,n2;
  /* drop off spaces */
  while (a<b && s[a]==' ') a++;
  /* skip numbers */
  i = a;
  while (i<b && s[i]!='-') i++;
  /* range between two numbers */
  if (i<b && s[i]=='-') {
    n1 = str_parse_unum_one(s,a,i+1);
    n2 = str_parse_unum_one(s,i+1,b);
    for (j=n1; j<=n2; j++)
      queue_uadd(q,j);
    }
  /* single number */
  else 
    queue_uadd(q,str_parse_unum_one(s,a,b));
  }

/* parse string containing list of numbers separated by commas and dashes 
 
   str - the string representing numbers
   num - number of numbers represented of the string (their values are rerurn
         in vector format) */
unsigned *str_parse_unum(char *str, unsigned *num) {
  unsigned i=0,start=0,end=str_length(str),*v=NULL;
  struct queue *q=queue_alloc();
  /* read numbers from string */
  while (i<end) {
    while (i<end && str[i]!=',')
      i++;
    str_parse_unum_commas(str,start,i,q);
    start = i+1;
    i++;
    }
  /* convert queue to vector */
  v = vec_ualloc(q->num);
  (*num) = q->num;
  for (i=0; i<(*num); i++)
    queue_uget(q,&(v[i]));
  /* finish */
  queue_free(q);
  return(v);
  }

/* -------------------------------------------------------------------------- */

/* read n-th word from the string and return pointer to it
 
   s - line with words
   c - word sequence number */
char *str_parse_word(char *s, unsigned c) {
  unsigned i,iw,n;
  if (!c)
    msg_error("column has to be specified by positive integer",1);
  i = 0;
  n = str_length(s);
  for (iw=0; iw<c; iw++) {
    while (i<n && (s[i]==' ' || s[i]=='\t')) i++;
    if ((iw+1)<c)
      while (i<n && (s[i]!=' ' && s[i]!='\t')) i++;
    if (i>=n)
      break;
    }
  if (iw<c)
    msg_error("specified column not found in the string",1);
  return(s+i);
  }

/* -------------------------------------------------------------------------- */
