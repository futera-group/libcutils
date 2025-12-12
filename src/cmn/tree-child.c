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
#include "cmn/tree.h"
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* increase child vector in tree graph node 
 
   p - pointer to the node
   n - pointer to new node */
void tree_child_vec_inc(struct tdata *p, struct tdata *n) {
  unsigned i;
  struct tdata **v;
  v = vec_talloc(p->t_num+1,sizeof(struct tdata*));
  for (i=0; i<p->t_num; i++)
    v[i] = p->t_child[i];
  vec_tfree(p->t_child);
  v[p->t_num] = n;
  p->t_child = v;
  p->t_num++;
  }

/* delete one child pointer from tree graph node 
 
   p - pointer to the node
   n - pointer to selected child */
void tree_child_vec_del(struct tdata *p, struct tdata *n) {
  unsigned i,j;
  struct tdata **v;
  if (!p->t_num)
    return;
  if (p->t_num==1) {
    if (p->t_child[0]==n) {
      p->t_child = vec_tfree(p->t_child);
      p->t_num--;
      }
    }
  else {
    for (i=0; i<p->t_num; i++)
      if (p->t_child[i]==n)
        break;
    if (i<p->t_num) {
      v = vec_talloc(p->t_num-1,sizeof(struct tdata*));
      for (j=0; j<i; j++)
        v[j] = p->t_child[j];
      for (j=i+1; j<p->t_num; j++)
        v[j-1] = p->t_child[j];
      p->t_child = vec_tfree(p->t_child);
      p->t_child = v;
      p->t_num--;
      }
    }
  }

/* add child node to parent one in tree graph
 
   t - pointer to the tree
   x - pointer to parent node
   r - type of data
   d - pointer to data */
struct tdata *tree_child_add(struct tree *t, struct tdata *x,
  short r, void *d) {
  struct tdata *p = NULL;
  /* allocate child */
  p = (struct tdata*)malloc(sizeof(struct tdata));
  if (!p)
    msg_error("cannot allocate child node of tree graph",1);
  p->t_child = NULL;
  p->t_parent = x;
  p->t_num = 0;
  p->t_id = 0;
  p->t_data = tree_alloc_data(r,d);
  /* connect to parent */
  tree_child_vec_inc(x,p);
  t->num++;
  return(p);
  }

/* sort childern in the tree node 
 
   p - the parent tree node
   f - comparing function */
void tree_child_sort(struct tdata *p, short(f)(void*,void*)) {
  short change;
  unsigned i,j;
  struct tdata *d;
  if (p->t_num>1) {
    /* bubble sort */
    for (i=1; i<p->t_num; i++) {
      change = 0;
      for (j=(p->t_num-1); j>=i; j--) {
        /* neighbour comparison */
        if (f(p->t_child[j-1]->t_data,p->t_child[j]->t_data)>0) {
          /* swap nodes */
          d = p->t_child[j-1];
          p->t_child[j-1] = p->t_child[j];
          p->t_child[j] = d;
          change++;
          }
        }
      /* sorted */
      if (!change)
        break;
      }
    }
  }

/* -------------------------------------------------------------------------- */
