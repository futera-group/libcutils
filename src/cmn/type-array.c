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
#include "cmn/message.h"
#include "cmn/types.h"

/* -------------------------------------------------------------------------- */

/* put one value into array of specified type
 
   v - pointer to the array
   n - ID of array element
   d - value for saving 
   t - type of data array */
void type_array_put(void *v, unsigned n, void *d, short t) {
  switch (t) {
    case TYPE_INT: 
      (*((int**)v))[n] = (*((int*)d)); break;
    case TYPE_UINT:
      (*((unsigned**)v))[n] = (*((unsigned*)d)); break;
    case TYPE_SINT:
      (*((short**)v))[n] = (*((short*)d)); break;
    case TYPE_USINT: 
      (*((unsigned short**)v))[n] = (*((unsigned short*)d)); break;
    case TYPE_LINT:
      (*((long**)v))[n] = (*((long*)d)); break;
    case TYPE_ULINT:
      (*((unsigned long**)v))[n] = (*((unsigned long*)d)); break;
    case TYPE_DOUBLE:
      (*((double**)v))[n] = (*((double*)d)); break;
    case TYPE_LDOUBLE:
      (*((long double**)v))[n] = (*((long double*)d)); break;
    case TYPE_STRING:
      msg_error("string arrays not supported in type_array_put",1);
    default:
      msg_error("unsupported type request passed to type_array_put",1);
      break;
    }
  }

/* -------------------------------------------------------------------------- */
