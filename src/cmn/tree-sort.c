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

/* sort descendants in subtree according to pathway lenghts
 
   t - pointer to the tree */
void tree_sort_len_sub(struct tdata *p) {
  unsigned i;
  for (i=0; i<p->t_num; i++)
    tree_sort_len_sub(p->t_child[i]);
  for (i=0; i<p->t_num; i++)
    if (p->t_child[i]->t_id>p->t_id)
      p->t_id = p->t_child[i]->t_id;
  }

/* sort descendants in tree according to pathway lenghts
 
   t - pointer to the tree */
void tree_sort_len(struct tree *t) {
  short check;
  unsigned i;
  struct tdata *p,*r;
  struct stack *s;
  /* set depth labels */
  tree_path_len_l(t,&p);
  /* reset labels to subtree lengths */
  tree_sort_len_sub(t->root);
  /* sort descendants */
  s = stack_alloc();
  stack_add(s,t->root);
  while (s->num) {
    p = stack_get(s);
    check = 1;
    while (check) {
      check = 0;
      for (i=1; i<p->t_num; i++)
        if (p->t_child[i-1]->t_id<p->t_child[i]->t_id) {
          r = p->t_child[i-1];
          p->t_child[i-1] = p->t_child[i];
          p->t_child[i] = r;
          check = 1;
          }
      for (i=0; i<p->t_num; i++)
        if (p->t_child[i]->t_num>1)
          stack_add(s,p->t_child[i]);
      }
    }
  stack_free(s);
  }

/* -------------------------------------------------------------------------- */
