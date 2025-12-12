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
#include "cmn/stack.h"
#include "cmn/string.h"
#include "cmn/tree.h"
#include "cmn/types.h"

/* -------------------------------------------------------------------------- */

/* allocate new tree structure */
struct tree *tree_alloc(void) {
  struct tree *t;
  t = (struct tree*)malloc(sizeof(struct tree));
  if (!t)
    msg_error("cannot allocate new tree",1);
  t->root = NULL;
  t->num = 0;
  return(t);
  }

/* allocate memory for tree node data and copy it
 
   t - type of data
   d - pointer to data */
void *tree_alloc_data(short t, void *d) {
  void *p = NULL;
  switch (t) {
    case TYPE_VOID:
      p = d;
      break;
    case TYPE_UINT:
      p = (unsigned*)malloc(sizeof(unsigned));
      if (!p)
        msg_error("cannot allocate memory for tree u-int data",1);
      (*((unsigned*)p))=(*((unsigned*)d));
      break;
    case TYPE_STRING:
      p = str_copy_new(d);
      break;
    default:
      msg_error("unsupported data type in tree_alloc_data",1);
    }
  return(p);
  }

/* clean tree structure

   t - pointer to the tree */
void tree_clean(struct tree *t) {
  unsigned i;
  struct tdata *p;
  struct stack *s;
  if (!t || !t->root)
    return;
  s = stack_alloc();
  stack_add(s,t->root);
  while (s->num) {
    p = stack_get(s);
    for (i=0; i<p->t_num; i++)
      stack_add(s,p->t_child[i]);
    free(p);
    }
  stack_free(s);
  t->root = NULL;
  t->num = 0;
  }

/* destroy the allocated tree

   t - pointer to the tree */
void tree_free(struct tree *t) {
  if (t) {
    tree_clean(t);
    free(t);
    }
  }

/* -------------------------------------------------------------------------- */
