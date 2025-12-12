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
#include "cmn/string.h"
#include "cmn/types.h"

/* -------------------------------------------------------------------------- */

/* allocate and copy value of of a given data type
 
   data - data which should be saved
   type - type of the data */
void *type_copy(void *data, short type) {
  unsigned n;
  void *d = NULL;
  switch (type) {
    case TYPE_INT: 
      d = malloc(sizeof(int));
      if (!d) msg_error("cannot allocate memory for list data (int)",1);
      (*((int*)d)) = (*((int*)data));
      break;
    case TYPE_UINT: 
      d = malloc(sizeof(unsigned));
      if (!d) msg_error("cannot allocate memory for list data (u-int)",1);
      (*((unsigned*)d)) = (*((unsigned*)data));
      break;
    case TYPE_SINT: 
      d = malloc(sizeof(short));
      if (!d) msg_error("cannot allocate memory for list data (s-int)",1);
      (*((short*)d)) = (*((short*)data));
      break;
    case TYPE_USINT: 
      d = malloc(sizeof(unsigned short));
      if (!d) msg_error("cannot allocate memory for list data (us-int)",1);
      (*((unsigned short*)d)) = (*((unsigned short*)data));
      break;
    case TYPE_LINT: 
      d = malloc(sizeof(long));
      if (!d) msg_error("cannot allocate memory for list data (l-int)",1);
      (*((long*)d)) = (*((long*)data));
      break;
    case TYPE_ULINT: 
      d = malloc(sizeof(unsigned long));
      if (!d) msg_error("cannot allocate memory for list data (ul-int)",1);
      (*((unsigned long*)d)) = (*((unsigned long*)data));
      break;
    case TYPE_DOUBLE: 
      d = malloc(sizeof(double));
      if (!d) msg_error("cannot allocate memory for list data (double)",1);
      (*((double*)d)) = (*((double*)data));
      break;
    case TYPE_LDOUBLE: 
      d = malloc(sizeof(long double));
      if (!d) msg_error("cannot allocate memory for list data (long double)",1);
      (*((long double*)d)) = (*((long double*)data));
      break;
    case TYPE_STRING: 
      n = str_length((char*)data);
      d = malloc((n+1)*sizeof(char));
      if (!d) msg_error("cannot allocate memory for list data (string)",1);
      str_copy((char*)d,(char*)data);
      ((char*)d)[n]='\0';
      break;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
