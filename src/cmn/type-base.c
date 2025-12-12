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
#include <stdarg.h>
#include <stdlib.h>
#include "cmn/message.h"
#include "cmn/types.h"

/* -------------------------------------------------------------------------- */

/* return string name of the specified data type

   t - the data type */
char* type_name(short t) {
  static char name[20] = "\0";
  switch (t) {
    case TYPE_VOID:    sprintf(name,"void");           break;
    case TYPE_INT:     sprintf(name,"int");            break;
    case TYPE_UINT:    sprintf(name,"unsigned");       break;
    case TYPE_SINT:    sprintf(name,"short");          break;
    case TYPE_USINT:   sprintf(name,"unsigned short"); break;
    case TYPE_LINT:    sprintf(name,"long");           break;
    case TYPE_ULINT:   sprintf(name,"unsigned long");  break;
    case TYPE_DOUBLE:  sprintf(name,"double");         break;
    case TYPE_LDOUBLE: sprintf(name,"long double");    break;
    case TYPE_STRING:  sprintf(name,"string");         break;
    default:           sprintf(name,"unknown");        break;
    }
  return(name);
  }

/* return string flag of the specified data type

   t - the data type */
char* type_flag(short t) {
  static char flag[20] = "\0";
  switch (t) {
    case TYPE_VOID:    sprintf(flag,"void");    break;
    case TYPE_INT:     sprintf(flag,"int");     break;
    case TYPE_UINT:    sprintf(flag,"uint");    break;
    case TYPE_SINT:    sprintf(flag,"sint");    break;
    case TYPE_USINT:   sprintf(flag,"usint");   break;
    case TYPE_LINT:    sprintf(flag,"lint");    break;
    case TYPE_ULINT:   sprintf(flag,"ulint");   break;
    case TYPE_DOUBLE:  sprintf(flag,"double");  break;
    case TYPE_LDOUBLE: sprintf(flag,"ldouble"); break;
    case TYPE_STRING:  sprintf(flag,"string");  break;
    default:           sprintf(flag,"unknown"); break;
    }
  return(flag);
  }

/* -------------------------------------------------------------------------- */

/* allocate memory for one specified data type
 
   d - data pointer
   t - type of data array */
void type_alloc(void *d, short t) {
  switch (t) {
    case TYPE_INT:
      (*((int**)d)) = 
        (int*)malloc(sizeof(int));
      if (!(*((int**)d)))
        msg_error("cannot allocate memory for int",1);
      break;
    case TYPE_UINT:
      (*((unsigned**)d)) = 
        (unsigned*)malloc(sizeof(unsigned));
      if (!(*((unsigned**)d)))
        msg_error("cannot allocate memory for u-int",1);
      break;
    case TYPE_SINT:
      (*((short**)d)) = 
        (short*)malloc(sizeof(short));
      if (!(*((short**)d)))
        msg_error("cannot allocate memory for s-int",1);
      break;
    case TYPE_USINT:
      (*((unsigned short**)d)) = 
        (unsigned short*)malloc(sizeof(unsigned short));
      if (!(*((unsigned short**)d)))
        msg_error("cannot allocate memory for us-int",1);
      break;
    case TYPE_LINT:
      (*((long**)d)) = 
        (long*)malloc(sizeof(long));
      if (!(*((long**)d)))
        msg_error("cannot allocate memory for l-int",1);
      break;
    case TYPE_ULINT:
      (*((unsigned long**)d)) = 
        (unsigned long*)malloc(sizeof(unsigned long));
      if (!(*((unsigned long**)d)))
        msg_error("cannot allocate memory for ul-int",1);
      break;
    case TYPE_DOUBLE:
      (*((double**)d)) = 
        (double*)malloc(sizeof(double));
      if (!(*((double**)d)))
        msg_error("cannot allocate memory for double",1);
      break;
    case TYPE_LDOUBLE:
      (*((long double**)d)) = 
        (long double*)malloc(sizeof(long double));
      if (!(*((long double**)d)))
        msg_error("cannot allocate memory for long double",1);
      break;
    case TYPE_STRING:
      msg_error("string allocation not supported by type_alloc",1);
      break;
    default:
      msg_error("unsupported type request passed to type_alloc",1);
      break;
    }
  }

/* allocate memory for vector of specified data types
 
   d - data pointer
   t - type of data array
   n - number of values */
void type_alloc_v(void *d, short t, unsigned n, ...) {
  unsigned i,m;
  va_list ap;
  va_start(ap,n);
  switch (t) {
    case TYPE_INT:
      (*((int**)d)) =
        (int*)malloc(n*sizeof(int));
      if (!(*((int**)d)))
        msg_error("cannot allocate memory for vector of int",1);
      break;
    case TYPE_UINT:
      (*((unsigned**)d)) =
        (unsigned*)malloc(n*sizeof(unsigned));
      if (!(*((unsigned**)d)))
        msg_error("cannot allocate memory for vector of u-int",1);
    case TYPE_SINT:
      (*((short**)d)) =
        (short*)malloc(n*sizeof(short));
      if (!(*((short**)d)))
        msg_error("cannot allocate memory for vector of s-int",1);
      break;
    case TYPE_USINT:
      (*((unsigned short**)d)) =
        (unsigned short*)malloc(n*sizeof(unsigned short));
      if (!(*((unsigned short**)d)))
        msg_error("cannot allocate memory for vector of us-int",1);
      break;
    case TYPE_LINT:
      (*((long**)d)) =
        (long*)malloc(n*sizeof(long));
      if (!(*((long**)d)))
        msg_error("cannot allocate memory for vector of l-int",1);
      break;
    case TYPE_ULINT:
      (*((unsigned long**)d)) =
        (unsigned long*)malloc(n*sizeof(unsigned long));
      if (!(*((unsigned long**)d)))
        msg_error("cannot allocate memory for vector of ul-int",1);
      break;
    case TYPE_DOUBLE:
      (*((double**)d)) =
        (double*)malloc(n*sizeof(double));
      if (!(*((double**)d)))
        msg_error("cannot allocate memory for vector of double",1);
      break;
    case TYPE_LDOUBLE:
      (*((long double**)d)) =
        (long double*)malloc(n*sizeof(long double));
      if (!(*((long double**)d)))
        msg_error("cannot allocate memory for vector of long double",1);
      break;
    case TYPE_STRING:
      m = va_arg(ap,unsigned);
      (*((char***)d)) =
        (char**)malloc(n*sizeof(char*));
      if (!(*((char***)d)))
        msg_error("cannot allocate memory for vector of strings",1);
      for (i=0; i<n; i++) {
        (*((char***)d))[i] = 
          (char*)malloc(m*sizeof(char));
        if (!((*((char***)d))[i]))
          msg_error("cannot allocate memory for vector of strings",1);
        }
      break;
    default:
      msg_error("unsupported type request passed to type_alloc_v",1);
      break;
    }
  va_end(ap);
  }

/* free memory allocated by type_alloc functions
 
   v - pointer to data storage */
void* type_free(void *v) {
  if (v)
    free(v);
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
