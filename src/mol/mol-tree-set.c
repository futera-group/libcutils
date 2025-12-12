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
#include <cmn/matrix.h>
#include <cmn/queue.h>
#include <cmn/stack.h>
#include <cmn/vector.h>
#include <cmn/tree.h>
#include <cmn/types.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* set tree marks for all atom in tree graph representation of molecule
 
   t - pointer to tree graph representation of molecule */
void mol_tree_set_marks(struct tree *t) {
  unsigned id;
  struct tdata *p,*r;
  struct stack *q;
  q = stack_alloc();
  stack_add(q,t->root);
  while (q->num) {
    p = stack_get(q);
    if (!p->t_parent)
      ((struct mol_tree_atom*)p->t_data)->tree = TREE_MARK_M;
    else if (!p->t_num)
      ((struct mol_tree_atom*)p->t_data)->tree = TREE_MARK_E;
    else {
      r = p->t_parent;
      for (id=0; id<r->t_num; id++)
        if (r->t_child[id]==p)
          break;
      if (!id && ((struct mol_tree_atom*)r->t_data)->tree==TREE_MARK_M)
        ((struct mol_tree_atom*)p->t_data)->tree = TREE_MARK_M;
      else {
        switch (p->t_num) {
          case 1: ((struct mol_tree_atom*)p->t_data)->tree = TREE_MARK_S; break;
          case 2: ((struct mol_tree_atom*)p->t_data)->tree = TREE_MARK_B; break;
          case 3: ((struct mol_tree_atom*)p->t_data)->tree = TREE_MARK_3; break;
          case 4: ((struct mol_tree_atom*)p->t_data)->tree = TREE_MARK_4; break;
          case 5: ((struct mol_tree_atom*)p->t_data)->tree = TREE_MARK_5; break;
          case 6: ((struct mol_tree_atom*)p->t_data)->tree = TREE_MARK_6; break;
          }
        }
      }
    for (id=0; id<p->t_num; id++)
      stack_add(q,p->t_child[id]);
    }
  stack_free(q);
  }

/* convert bond molecular representation to graph tree
 
   b  - bond array with atom IDs
   nb - number of bonds
   na - number of atoms
   l  - pointer to loop array (output)
   nl - number of loops */
struct tree *mol_tree_set(unsigned **b, unsigned nb, unsigned na,
  unsigned ***l, unsigned *nl) {
  unsigned i,rt,*lp,id = 0;
  short *at,*bn;
  struct mol_tree_atom *x;
  struct tdata *p,*d = NULL;
  struct queue *q,*o;
  struct tree *t;
  /* initialization */
  t = tree_alloc();
  at = vec_sialloc(na);
  vec_siset(at,0,na);
  bn = vec_sialloc(nb);
  vec_siset(bn,0,nb);
  q = queue_alloc();
  o = queue_alloc();
  /* set root */
  x = mol_tree_atom_new();
  x->id = b[0][0];
  p = tree_root_set(t,TYPE_VOID,x);
  queue_add(q,p);
  at[b[0][0]] = 1;
  x = mol_tree_atom_new();
  x->id = b[0][1];
  p = tree_child_add(t,p,TYPE_VOID,x);
  queue_add(q,p);
  at[b[0][1]] = 1;
  bn[0] = 1;
  /* connect all other atoms */
  while (q->num) {
    p = queue_get(q);
    rt = ((struct mol_tree_atom*)p->t_data)->id;
    for (i=0; i<nb; i++)
      if (!bn[i] && (b[i][0]==rt || b[i][1]==rt)) {
        if (at[b[i][0]] && at[b[i][1]]) {
          lp = vec_ualloc(2);
          lp[0] = b[i][0];
          lp[1] = b[i][1];
          queue_add(o,lp);
          }
        else {
          x = mol_tree_atom_new();
          x->id = (b[i][0]==rt ? b[i][1] : b[i][0]);
          d = tree_child_add(t,p,TYPE_VOID,x);
          at[x->id] = 1;
          queue_add(q,d);
          }
        bn[i] = 1;
        }
    }
  /* create loop array */
  (*nl) = o->num;
  (*l) = mat_ualloc(o->num,2);
  while (o->num) {
    lp = queue_get(o);
    (*l)[id][0] = lp[0];
    (*l)[id][1] = lp[1];
    free(lp);
    id++;
    }
  /* clean memory */
  vec_sifree(at);
  vec_sifree(bn);
  queue_free(q);
  queue_free(o);
  /* sort children */
  tree_root_maxpath(t);
  tree_sort_len(t);
  /* set tree marks */
  mol_tree_set_marks(t);
  return(t);
  }

/* -------------------------------------------------------------------------- */
