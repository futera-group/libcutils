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
#include "cmn/stack.h"
#include "cmn/tree.h"

/* -------------------------------------------------------------------------- */

/* print out string data saved in the tree structure

   t - pointer to the tree */
void tree_sprint(struct tree *t) {
  unsigned i;
  struct stack *q;
  struct tdata *p;
  if (t) {
    /* initialization */
    q = stack_alloc();
    if (t->root) {
      t->root->t_id = 0;
      stack_add(q,t->root);
      }
    /* print tree */
    while (q->num) {
      p = stack_get(q);
      if (p->t_data) {
        if (p->t_id)
          printf("%*c%s\n",2*p->t_id,' ',(char*)p->t_data);
        else
          printf("%s\n",(char*)p->t_data);
        }
      for (i=0; i<p->t_num; i++) {
        p->t_child[p->t_num-i-1]->t_id = p->t_id+1;
        stack_add(q,p->t_child[p->t_num-i-1]);
        }
      }
    stack_free(q);
    }
  }

/* -------------------------------------------------------------------------- */
