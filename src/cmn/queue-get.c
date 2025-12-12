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
#include <string.h>
#include "cmn/queue.h"

/* -------------------------------------------------------------------------- */

/* return pointer to the first data item in the queue 
 
   q - pointer to the queue */
void *queue_get(struct queue *q) {
  void *data=NULL;
  struct qdata *qp;
  if (q->num>0) {
    qp = q->data;
    q->data = qp->q_next;
    data = qp->q_data;
    free(qp);
    (q->num)--;
    }
  return(data);
  }

/* -------------------------------------------------------------------------- */

/* get integer from the queue
  
   q   - pointer to the queue
   val - value which will be get from the queue */
long unsigned queue_iget(struct queue *q, int *val) {
  void *data;
  data = queue_get(q);
  (*val) = *(int*)data;
  free(data);
  return(q->num);
  }

/* get unsigned int from the queue
  
   q   - pointer to the queue
   val - value which will be get from the queue */
long unsigned queue_uget(struct queue *q, unsigned *val) {
  void *data;
  data = queue_get(q);
  (*val) = *(unsigned*)data;
  free(data);
  return(q->num);
  }

/* get short integer from the queue
  
   q   - pointer to the queue
   val - value which will be get from the queue */
long unsigned queue_siget(struct queue *q, short *val) {
  void *data;
  data = queue_get(q);
  (*val) = *(short*)data;
  free(data);
  return(q->num);
  }

/* get short unsigned integer from the queue
  
   q   - pointer to the queue
   val - value which will be get from the queue */
long unsigned queue_suget(struct queue *q, short unsigned *val) {
  void *data;
  data = queue_get(q);
  (*val) = *(short unsigned*)data;
  free(data);
  return(q->num);
  }

/* get long integer from the queue
  
   q   - pointer to the queue
   val - value which will be get from the queue */
long unsigned queue_liget(struct queue *q, long *val) {
  void *data;
  data = queue_get(q);
  (*val) = *(long*)data;
  free(data);
  return(q->num);
  }

/* get long unsigned integer from the queue
  
   q   - pointer to the queue
   val - value which will be get from the queue */
long unsigned queue_luget(struct queue *q, long unsigned *val) {
  void *data;
  data = queue_get(q);
  (*val) = *(long unsigned*)data;
  free(data);
  return(q->num);
  }

/* -------------------------------------------------------------------------- */

/* get double precision real from the queue
  
   q   - pointer to the queue
   val - value which will be get from the queue */
long unsigned queue_fget(struct queue *q, double *val) {
  void *data;
  data = queue_get(q);
  (*val) = *(double*)data;
  free(data);
  return(q->num);
  }

/* get long double precision real from the queue
  
   q   - pointer to the queue
   val - value which will be get from the queue */
long unsigned queue_lfget(struct queue *q, long double *val) {
  void *data;
  data = queue_get(q);
  (*val) = *(long double*)data;
  free(data);
  return(q->num);
  }

/* -------------------------------------------------------------------------- */

/* get integer from the queue
  
   q   - pointer to the queue
   val - value which will be get from the queue */
long unsigned queue_cget(struct queue *q, char *val) {
  void *data;
  data = queue_get(q);
  (*val) = *(char*)data;
  free(data);
  return(q->num);
  }

/* get string from the queue and return actual number of items in the queue
 
   q   - pointer to the queue
   str - the string */
long unsigned queue_sget(struct queue *q, char *str) {
  void *data;
  data = queue_get(q);
  strcpy(str,(char*)data);
  free(data);
  return(q->num);
  }

/* -------------------------------------------------------------------------- */
