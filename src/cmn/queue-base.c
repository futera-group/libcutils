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

/* -------------------------------------------------------------------------- */

/* allocate new queue structure */
struct queue *queue_alloc(void) {
  struct queue *q;
  q = (struct queue*)malloc(sizeof(struct queue));
  if (!q)
    msg_error("cannot allocate new queue",1);
  q->data = NULL;
  q->num = 0;
  return(q);
  }

/* clean queue structure
 
   q - pointer to the queue */
void queue_clean(struct queue *q) {
  struct qdata *qp,*qp2;
  if (!q) 
    return;
  qp = q->data;
  while (qp) {
    qp2 = qp;
    qp = qp->q_next;
    free(qp2->q_data);
    free(qp2);
    }
  q->num = 0;
  }

/* destroy the allocated queue
 
   q - pointer to the queue */
void queue_free(struct queue *q) {
  if (q) {
    queue_clean(q);
    free(q);
    }
  }

/* -------------------------------------------------------------------------- */
