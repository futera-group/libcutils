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

/* return pointer to specific element in a rank-1 data array

   d    - pointer to data value
   type - type data
   id   - ID of the array element
   n    - first dimension of the array */
void *cfg_rank1_pointer(void *d, short type, unsigned id, unsigned n) {
  void *p = NULL;
  switch (type) {
    case TYPE_INT: /* integer */
      p = (*((int**)d))+id;
      break;
    case TYPE_UINT: /* unsigned integer */
      p = (*((unsigned**)d))+id;
      break;
    case TYPE_SINT: /* short integer */
      p = (*((short**)d))+id;
      break;
    case TYPE_USINT: /* unsigned short integer */
      p = (*((unsigned short**)d))+id;
      break;
    case TYPE_LINT: /* long integer */
      p = (*((long**)d))+id;
      break;
    case TYPE_ULINT: /* unsigned long integer */
      p = (*((unsigned long**)d))+id;
      break;
    case TYPE_DOUBLE: /* double-precision real */
      p = (*((double**)d))+id;
      break;
    case TYPE_LDOUBLE: /* long double-precision real */
      p = (*((long double**)d))+id;
      break;
    case TYPE_STRING: /* string */
      p = (*((char***)d))+id;
      break;
    default:
      msg_error_f("attempt to get rank-1 pointer of unknown"
        " data type (%d)",1,type);
    }
  return(p);
  }

/* return pointer to specific element in a rank-2 data array

   d    - pointer to data value
   type - type data
   id   - ID of the array element
   n    - first dimension of the array */
void *cfg_rank2_pointer(void *d, short type, unsigned id, unsigned n) {
  void *p = NULL;
  switch (type) {
    case TYPE_INT: /* integer */
      p = (*((int***)d))+id;
      break;
    case TYPE_UINT: /* unsigned integer */
      p = (*((unsigned***)d))+id;
      break;
    case TYPE_SINT: /* short integer */
      p = (*((short***)d))+id;
      break;
    case TYPE_USINT: /* unsigned short integer */
      p = (*((unsigned short***)d))+id;
      break;
    case TYPE_LINT: /* long integer */
      p = (*((long***)d))+id;
      break;
    case TYPE_ULINT: /* unsigned long integer */
      p = (*((unsigned long***)d))+id;
      break;
    case TYPE_DOUBLE: /* double-precision real */
      p = (*((double***)d))+id;
      break;
    case TYPE_LDOUBLE: /* long double-precision real */
      p = (*((long double***)d))+id;
      break;
    case TYPE_STRING: /* string */
      p = (*((char****)d))+id;
      break;
    default:
      msg_error_f("attempt to get rank-2 pointer of unknown"
        " data type (%d)",1,type);
    }
  return(p);
  }

/* return pointer to specific element in a rank-3 data array

   d    - pointer to data value
   type - type data
   id   - ID of the array element
   n    - first dimension of the array */
void *cfg_rank3_pointer(void *d, short type, unsigned id, unsigned n) {
  void *p = NULL;
  switch (type) {
    case TYPE_INT: /* integer */
      p = (*((int****)d))+id;
      break;
    case TYPE_UINT: /* unsigned integer */
      p = (*((unsigned****)d))+id;
      break;
    case TYPE_SINT: /* short integer */
      p = (*((short****)d))+id;
      break;
    case TYPE_USINT: /* unsigned short integer */
      p = (*((unsigned short****)d))+id;
      break;
    case TYPE_LINT: /* long integer */
      p = (*((long****)d))+id;
      break;
    case TYPE_ULINT: /* unsigned long integer */
      p = (*((unsigned long****)d))+id;
      break;
    case TYPE_DOUBLE: /* double-precision real */
      p = (*((double****)d))+id;
      break;
    case TYPE_LDOUBLE: /* long double-precision real */
      p = (*((long double****)d))+id;
      break;
    case TYPE_STRING: /* string */
      p = (*((char*****)d))+id;
      break;
    default:
      msg_error_f("attempt to get rank-3 pointer of unknown"
        " data type (%d)",1,type);
    }
  return(p);
  }

/* return pointer to specific element in a rank-4 data array

   d    - pointer to data value
   type - type data
   id   - ID of the array element
   n    - first dimension of the array */
void *cfg_rank4_pointer(void *d, short type, unsigned id, unsigned n) {
  void *p = NULL;
  switch (type) {
    case TYPE_INT: /* integer */
      p = (*((int*****)d))+id;
      break;
    case TYPE_UINT: /* unsigned integer */
      p = (*((unsigned*****)d))+id;
      break;
    case TYPE_SINT: /* short integer */
      p = (*((short*****)d))+id;
      break;
    case TYPE_USINT: /* unsigned short integer */
      p = (*((unsigned short*****)d))+id;
      break;
    case TYPE_LINT: /* long integer */
      p = (*((long*****)d))+id;
      break;
    case TYPE_ULINT: /* unsigned long integer */
      p = (*((unsigned long*****)d))+id;
      break;
    case TYPE_DOUBLE: /* double-precision real */
      p = (*((double*****)d))+id;
      break;
    case TYPE_LDOUBLE: /* long double-precision real */
      p = (*((long double*****)d))+id;
      break;
    case TYPE_STRING: /* string */
      p = (*((char******)d))+id;
      break;
    default:
      msg_error_f("attempt to get rank-4 pointer of unknown"
        " data type (%d)",1,type);
    }
  return(p);
  }

/* return pointer to specific element in a data array
 
   d    - pointer to data value
   type - type data
   rank - rank of the data array
   id   - ID of the array element
   n    - first dimension of the array */
void *cfg_data_pointer(void *d, short type, int rank,
  unsigned id, unsigned n) {
  void *p = NULL;
  if (d && rank>0) {
    if (id>=n)
      msg_error_f("attempt to get pointer of invalid element (%d/%d)",1,id+1,n);
    switch (rank) {
      case 1: p = cfg_rank1_pointer(d,type,id,n); break;
      case 2: p = cfg_rank2_pointer(d,type,id,n); break;
      case 3: p = cfg_rank3_pointer(d,type,id,n); break;
      case 4: p = cfg_rank4_pointer(d,type,id,n); break;
      default:
        msg_error_f("attempt to get pointer of unsupported"
          " data rank (%d)",1,rank);
      }
    }
  return(p);
  }

/* -------------------------------------------------------------------------- */
