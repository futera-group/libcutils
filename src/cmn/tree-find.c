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
#include "cmn/stack.h"
#include "cmn/tree.h"
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* find element in tree according to specified comparing function
  
   t - the tree data structure
   k - data to search
   f - the comparing function
   n - number of found records (output) */
void *tree_find_data_t(struct tree *t, void *k,
  short(f)(void*,void*), unsigned *n) {
  unsigned i,id;
  struct tdata *p,**d = NULL;
  struct queue *q;
  struct stack *s;
  (*n) = 0;
  if (!t || !t->root)
    return(d);
  s = stack_alloc();
  q = queue_alloc();
  /* go through the tree */
  stack_add(s,t->root);
  while (s->num) {
    p = stack_get(s);
    for (i=0; i<p->t_num; i++)
      stack_add(s,p->t_child[i]);
    /* save results to queue */
    if (f(p->t_data,k))
      queue_add(q,p);
    }
  stack_free(s);
  /* convert queue to array */
  id = 0;
  (*n) = q->num;
  d = vec_talloc(sizeof(struct tdata*),q->num);
  while (q->num) {
    p = queue_get(q);
    d[id] = p;
    id++;
    }
  queue_free(q);
  return(d);
  }

/* -------------------------------------------------------------------------- */
