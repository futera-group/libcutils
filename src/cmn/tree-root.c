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
#include "cmn/tree.h"

/* -------------------------------------------------------------------------- */

/* set root of tree structure
 
   t - pointer to the tree
   r - type of data
   d - pointer to data */
struct tdata *tree_root_set(struct tree *t, short r, void *d) {
  struct tdata *p = NULL;
  p = (struct tdata*)malloc(sizeof(struct tdata));
  if (!p)
    msg_error("cannot allocate root node of tree graph",1);
  p->t_child = NULL;
  p->t_parent = NULL;
  p->t_num = 0;
  p->t_id = 0;
  p->t_data = tree_alloc_data(r,d);
  t->root = p;
  t->num = 1;
  return(p);
  }

/* change root of the tree to specifed node
 
   t - pointer to the tree
   r - pointer to the new root */
void tree_root_change(struct tree *t, struct tdata *r) {
  struct tdata *p = r,*f,*l = NULL;
  while (p->t_parent) {
    f = p->t_parent;
    tree_child_vec_inc(p,f);
    tree_child_vec_del(f,p);
    p->t_parent = l;
    l = p;
    p = f;
    }
  p->t_parent = l;
  t->root = r;
  }

/* reorganize tree to maximalize path from root to leaf
 
   t - pointer to the tree */
void tree_root_maxpath(struct tree *t) {
  unsigned i,ln,ln_l = 0;
  struct tdata *p,*d,*p_l = t->root;
  struct stack *s1,*s2;
  s1 = stack_alloc();
  s2 = stack_alloc();
  stack_add(s1,t->root);
  stack_add(s2,t->root);
  while (s1->num) {
    p = stack_get(s1);
    for (i=0; i<p->t_num; i++) {
      stack_add(s1,p->t_child[i]);
      stack_add(s2,p->t_child[i]);
      }
    }
  stack_free(s1);
  while (s2->num) {
    p = stack_get(s2);
    tree_root_change(t,p);
    ln=tree_path_len_l(t,&d);
    if (ln>ln_l) {
      ln_l = ln;
      p_l = p;
      }
    }
  stack_free(s2);
  tree_root_change(t,p_l);
  }

/* -------------------------------------------------------------------------- */
