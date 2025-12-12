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
#include <stdlib.h>
#include "cmn/message.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* return lenght of string representing int

  num = integer */
unsigned str_inum_len(int num) {
  int n;
  unsigned len;
  n = abs(num);
  len = 0;
  while (n>0) {
    n /= 10;
    len++;
    }
  if (num<0)
    len++;
  return(len);
  }

/* return lenght of string representing unsigned int

  num = unsigned int */
unsigned str_unum_len(unsigned num) {
  unsigned len,n;
  n = num;
  len = 0;
  while (n>0) {
    n /= 10;
    len++;
    }
  return(len);
  }

/* return lenght of string representing short int

  num = short integer */
unsigned str_sinum_len(short int num) {
  short int n;
  unsigned len;
  n = abs(num);
  len = 0;
  while (n>0) {
    n /= 10;
    len++;
    }
  if (num<0)
    len++;
  return(len);
  }

/* return lenght of string representing unsigned short

  num = unsigned short */
unsigned str_sunum_len(unsigned short num) {
  unsigned len;
  unsigned short n;
  n = num;
  len = 0;
  while (n>0) {
    n /= 10;
    len++;
    }
  return(len);
  }

/* return lenght of string representing long int

  num = long integer */
unsigned str_linum_len(long int num) {
  long int n;
  unsigned len;
  n = labs(num);
  len = 0;
  while (n>0) {
    n /= 10;
    len++;
    }
  if (num<0)
    len++;
  return(len);
  }

/* return lenght of string representing unsigned long

  num = unsigned long */
unsigned str_lunum_len(unsigned long num) {
  unsigned len;
  unsigned long n;
  n = num;
  len = 0;
  while (n>0) {
    n /= 10;
    len++;
    }
  return(len);
  }

/* -------------------------------------------------------------------------- */

/* convert vector of numbers to string representation (u-int)

   str - string where range of numbers is writen
   vec - vector of numbers
   num - dimension of the vector */
void str_unum_range(char *str, unsigned *vec, unsigned num) {
  unsigned i=0,n=0,first=0,last=0;
  while (i<num) {
    if (first>=num)
      break;
    for (i=first+1; i<num; i++) {
      if (vec[i]!=vec[i-1]+1) {
        last = i-1;
        break;
        }
      if (i==num-1) {
        last = i;
        break;
        }
      }
    sprintf(str+n,"%d",vec[first]);
    n += str_unum_len(vec[first]);
    if (last>first) {
      sprintf(str+n,"-%d",vec[last]);
      n += (1+str_unum_len(vec[last]));
      }
    if (last<(num-1)) {
      sprintf(str+n,",");
      n++;
      }
    first = last+1;
    last = last+1;
    }
  str[n]='\0';
  }

/* convert vector of c-indices to string representation (u-int)

   str - string where range of numbers is writen
   vec - vector of numbers
   num - dimension of the vector */
void str_unum_Range(char *str, unsigned *vec, unsigned num) {
  unsigned i;
  for (i=0; i<num; i++)
    vec[i]++;
  str_unum_range(str,vec,num);
  for (i=0; i<num; i++)
    vec[i]--;
  }

/* -------------------------------------------------------------------------- */

/* convert fortran string number to double precision real

   snum - string number if fortran format (D instead of E in exponents) */
double str_fnum_d2e(char *snum) {
  unsigned i;
  double dnum;
  if (str_length(snum)<2)
    return(0.0);
  for (i=0; i<str_length(snum); i++) {
    if (snum[i]=='d')
      snum[i]='e';
    else if (snum[i]=='D')
      snum[i]='E';
    }
  if (sscanf(snum,"%lf",&dnum)!=1)
    msg_error_f("cannot convert fortran double precision type \"%s\"",1,snum);
  return(dnum);
  }

/* -------------------------------------------------------------------------- */

/* convert unsigned integer number to string representation

   n - the number */
char* str_unum_i2s_new(unsigned n) {
  char *s = NULL;
  s = str_new(str_unum_len(n));
  sprintf(s,"%d",n);
  return(s);
  }

/* -------------------------------------------------------------------------- */
