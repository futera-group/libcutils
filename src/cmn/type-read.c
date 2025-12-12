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

#include <limits.h>
#include <stdlib.h>
#include "cmn/message.h"
#include "cmn/string.h"
#include "cmn/types.h"

/* -------------------------------------------------------------------------- */

/* read value of specified type from string 
 
   s - the string
   d - data pointer
   t - type of data array */
short type_read(char *s, void *d, short t) {
  unsigned long lu;
  long li;
  char *cp;
  cp = s;
  switch (t) {
    /* signed integers */
    case TYPE_INT:
    case TYPE_SINT:
    case TYPE_LINT:
      li = strtol(s,&cp,10);
      if (cp==s) 
        return(0);
      if (t==TYPE_SINT) {
        if (li>SHRT_MAX || li<SHRT_MIN) 
          return(0);
        (*(short*)d) = (short)li;
        }
      else if (t==TYPE_INT) {
        if (li>INT_MAX || li<INT_MIN)  
          return(0);
        (*(int*)d) = (int)li;
        }
      else
        (*(long*)d) = li;
      break;
    /* signed integers */
    case TYPE_UINT:
    case TYPE_USINT:
    case TYPE_ULINT:
      lu = strtoul(s,&cp,10);
      if (cp==s)
        return(0);
      if (t==TYPE_USINT) {
        if (lu>USHRT_MAX)
          return(0);
        (*(unsigned short*)d) = (unsigned short)lu;
        }
      else if (t==TYPE_UINT) {
        if (lu>UINT_MAX)
          return(0);
        (*(unsigned*)d) = (unsigned)lu;
        }
      else
        (*(unsigned long*)d) = lu;
      break;
    /* double-precision real */
    case TYPE_DOUBLE:
      (*(double*)d) = strtod(s,&cp);
      if (cp==s)
        return(0);
      break;
    case TYPE_LDOUBLE:
      (*(long double*)d) = strtold(s,&cp);
      if (cp==s)
        return(0);
      break;
    /* string */
    case TYPE_STRING:
      (*((char**)d)) = str_copy_new(s);
      break;
    default:
      msg_error("unsupported type request passed to type_read",1);
      break;
    }
  return(1);
  }

/* -------------------------------------------------------------------------- */
