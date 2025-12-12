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
#include "cmn/list.h"
#include "cmn/types.h"
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* convert data saved in a list to vector
 
   q - the list data 
   t - data type
   n - number of converted items */
void *list_conv_vec(struct list *q, short t, long unsigned *n) {
  long unsigned i=0;
  void *v = NULL;
  (*n) = q->num;
  switch (t) {
    case TYPE_INT:
      v = vec_ialloc(q->num);
      while (q->num) {
        ((int*)v)[i++] = *((int*)list_get_begin(q));
        list_del_begin(q,free);
        }
      break;
    case TYPE_UINT:
      v = vec_ualloc(q->num);
      while (q->num) {
        ((unsigned*)v)[i++] = *((unsigned*)list_get_begin(q));
        list_del_begin(q,free);
        }
      break;
    case TYPE_SINT:
      v = vec_sialloc(q->num);
      while (q->num) {
        ((short*)v)[i++] = *((short*)list_get_begin(q));
        list_del_begin(q,free);
        }
      break;
    case TYPE_USINT:
      v = vec_sualloc(q->num);
      while (q->num) {
        ((short unsigned*)v)[i++] = *((short unsigned*)list_get_begin(q));
        list_del_begin(q,free);
        }
      break;
    case TYPE_LINT:
      v = vec_lialloc(q->num);
      while (q->num) {
        ((long*)v)[i++] = *((long*)list_get_begin(q));
        list_del_begin(q,free);
        }
      break;
    case TYPE_ULINT: 
      v = vec_lualloc(q->num);
      while (q->num) {
        ((unsigned long*)v)[i++] = *((unsigned long*)list_get_begin(q));
        list_del_begin(q,free);
        }
      break;
    case TYPE_DOUBLE:
      v = vec_falloc(q->num);
      while (q->num) {
        ((double*)v)[i++] = *((double*)list_get_begin(q));
        list_del_begin(q,free);
        }
      break;
    case TYPE_LDOUBLE:
      v = vec_lfalloc(q->num);
      while (q->num) {
        ((long double*)v)[i++] = *((long double*)list_get_begin(q));
        list_del_begin(q,free);
        }
      break;
    }
  return(v);
  }

/* -------------------------------------------------------------------------- */
