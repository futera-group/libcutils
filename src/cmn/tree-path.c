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

#include "cmn/stack.h"
#include "cmn/tree.h"

/* -------------------------------------------------------------------------- */

/* return length of path from root to given node
 
   t - pointer to the tree
   f - pointer to the node */
unsigned tree_path_len_p(struct tree *t, struct tdata *f) {
  unsigned len=0;
  struct tdata *p;
  p = f;
  while (p->t_parent) {
    p = p->t_parent;
    len++;
    }
  return(len);
  }

/* return length of longest path in the tree graph
 
   t - pointer to the tree
   f - pointer to most distant leaf (output) */
unsigned tree_path_len_l(struct tree *t, struct tdata **f) {
  unsigned i,len=0;
  struct stack *s;
  struct tdata *p;
  if (f)
    (*f) = t->root;
  s = stack_alloc();
  stack_add(s,t->root);
  while (s->num) {
    p = stack_get(s);
    if (p->t_parent)
      p->t_id = (p->t_parent->t_id+1);
    else 
      p->t_id = 0;
    for (i=0; i<p->t_num; i++)
      stack_add(s,p->t_child[i]);
    if (!p->t_num && p->t_id>len) {
      len = p->t_id;
      if (f)
        (*f) = p;
      }
    }
  stack_free(s);
  return(len);
  }

/* return length of shortest path in the tree graph
 
   t - pointer to the tree
   f - pointer to most distant leaf (output) */
unsigned tree_path_len_s(struct tree *t, struct tdata **f) {
  unsigned i,len=t->num;
  struct stack *s;
  struct tdata *p;
  if (f)
    (*f) = t->root;
  s = stack_alloc();
  stack_add(s,t->root);
  while (s->num) {
    p = stack_get(s);
    if (p->t_parent)
      p->t_id = (p->t_parent->t_id+1);
    else 
      p->t_id = 0;
    for (i=0; i<p->t_num; i++)
      stack_add(s,p->t_child[i]);
    if (!p->t_num && p->t_id<len) {
      len = p->t_id;
      if (f)
        (*f) = p;
      }
    }
  stack_free(s);
  return(len);
  }

/* -------------------------------------------------------------------------- */
