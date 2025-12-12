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
#include "cmn/queue.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* add new data to the queue
   
   q    - pointer to the queue
   data - pointer to the data */
void queue_add(struct queue *q, void *data) {
  struct qdata *qp=NULL,*qtp;
  qp = malloc(sizeof(struct qdata));
  if (!qp)
    msg_error("cannot allocate memory for queue cell",1);
  qp->q_data = data;
  qp->q_next = NULL;
  if (!q->num)
    q->data=qp;
  else {
    qtp = q->data;
    while (qtp->q_next)
      qtp = qtp->q_next;
    qtp->q_next = qp;
    }
  (q->num)++;
  }

/* -------------------------------------------------------------------------- */

/* add integer to the queue
  
   q   - pointer to the queue
   val - value which will be inserted to the queue */
long unsigned queue_iadd(struct queue *q, int val) {
  int *data;
  data = (int*)malloc(sizeof(int));
  if (!data)
    msg_error("cannot allocate memory for i-queue data",1);
  (*data) = val;
  queue_add(q,data);
  return(q->num);
  }

/* add unsigned int to the queue
  
   q   - pointer to the queue
   val - value which will be inserted to the queue */
long unsigned queue_uadd(struct queue *q, unsigned val) {
  unsigned *data;
  data = (unsigned*)malloc(sizeof(unsigned));
  if (!data)
    msg_error("cannot allocate memory for u-queue data",1);
  (*data) = val;
  queue_add(q,data);
  return(q->num);
  }

/* add short integer to the queue
  
   q   - pointer to the queue
   val - value which will be inserted to the queue */
long unsigned queue_siadd(struct queue *q, short val) {
  short *data;
  data = (short*)malloc(sizeof(short));
  if (!data)
    msg_error("cannot allocate memory for si-queue data",1);
  (*data) = val;
  queue_add(q,data);
  return(q->num);
  }

/* add short unsigned integer to the queue
  
   q   - pointer to the queue
   val - value which will be inserted to the queue */
long unsigned queue_suadd(struct queue *q, short unsigned val) {
  short unsigned *data;
  data = (short unsigned*)malloc(sizeof(short unsigned));
  if (!data)
    msg_error("cannot allocate memory for su-queue data",1);
  (*data) = val;
  queue_add(q,data);
  return(q->num);
  }

/* add long integer to the queue
  
   q   - pointer to the queue
   val - value which will be inserted to the queue */
long unsigned queue_liadd(struct queue *q, long val) {
  long *data;
  data = (long*)malloc(sizeof(long));
  if (!data)
    msg_error("cannot allocate memory for li-queue data",1);
  (*data) = val;
  queue_add(q,data);
  return(q->num);
  }

/* add long unsigned integer to the queue
  
   q   - pointer to the queue
   val - value which will be inserted to the queue */
long unsigned queue_luadd(struct queue *q, long unsigned val) {
  long unsigned *data;
  data = (long unsigned*)malloc(sizeof(long unsigned));
  if (!data)
    msg_error("cannot allocate memory for lu-queue data",1);
  (*data) = val;
  queue_add(q,data);
  return(q->num);
  }

/* -------------------------------------------------------------------------- */

/* add double precision real to the queue
  
   q   - pointer to the queue
   val - value which will be inserted to the queue */
long unsigned queue_fadd(struct queue *q, double val) {
  double *data;
  data = (double*)malloc(sizeof(double));
  if (!data)
    msg_error("cannot allocate memory for f-queue data",1);
  (*data) = val;
  queue_add(q,data);
  return(q->num);
  }

/* add long double precision real to the queue
  
   q   - pointer to the queue
   val - value which will be inserted to the queue */
long unsigned queue_lfadd(struct queue *q, long double val) {
  long double *data;
  data = (long double*)malloc(sizeof(long double));
  if (!data)
    msg_error("cannot allocate memory for Lf-queue data",1);
  (*data) = val;
  queue_add(q,data);
  return(q->num);
  }

/* -------------------------------------------------------------------------- */

/* add character to the queue
  
   q   - pointer to the queue
   val - value which will be inserted to the queue */
long unsigned queue_cadd(struct queue *q, char val) {
  char *data;
  data = (char*)malloc(sizeof(char));
  if (!data)
    msg_error("cannot allocate memory for c-queue data",1);
  (*data) = val;
  queue_add(q,data);
  return(q->num);
  }

/* add string to the queue and return actual number of items in the queue
   
   q   - pointer to the queue
   str - the string */
long unsigned queue_sadd(struct queue *q, char *str) {
  queue_add(q,str_copy_new(str));
  return(q->num);
  }


/* -------------------------------------------------------------------------- */
