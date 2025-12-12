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
#include "cmn/queue.h"
#include "cmn/types.h"
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* convert data saved in a queue to vector
 
   q - the queue data 
   t - data type
   n - number of converted items */
void *queue_conv_vec(struct queue *q, short t, long unsigned *n) {
  long unsigned i=0;
  void *v = NULL;
  (*n) = q->num;
  switch (t) {
    case TYPE_INT:
      v = vec_ialloc(q->num);
      while (q->num)
        queue_iget(q,&(((int*)v)[i++]));
      break;
    case TYPE_UINT:
      v = vec_ualloc(q->num);
      while (q->num)
        queue_uget(q,&(((unsigned*)v)[i++]));
      break;
    case TYPE_SINT:
      v = vec_sialloc(q->num);
      while (q->num)
        queue_siget(q,&(((short*)v)[i++]));
      break;
    case TYPE_USINT:
      v = vec_sualloc(q->num);
      while (q->num)
        queue_suget(q,&(((short unsigned*)v)[i++]));
      break;
    case TYPE_LINT:
      v = vec_lialloc(q->num);
      while (q->num)
        queue_liget(q,&(((long*)v)[i++]));
      break;
    case TYPE_ULINT: 
      v = vec_lualloc(q->num);
      while (q->num)
        queue_luget(q,&(((unsigned long*)v)[i++]));
      break;
    case TYPE_DOUBLE:
      v = vec_falloc(q->num);
      while (q->num)
        queue_fget(q,&(((double*)v)[i++]));
      break;
    }
  return(v);
  }

/* -------------------------------------------------------------------------- */
