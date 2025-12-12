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
#include "cmn/config.h"
#include "cmn/list.h"
#include "cmn/queue.h"
#include "cmn/string.h"
#include "cmn/types.h"

/* -------------------------------------------------------------------------- */

/* define tree of sub-item templates in the new generated suboption

   d  - main config data struct
   p  - parent template item
   c  - the new generated suboption
   id - replacement ID of the first dependency mark
   n  - first dimension of the first dependency array */
void cfg_set_subsubitems(struct cfg *d, struct cfg_item *p,
  struct cfg_item *c, unsigned id, unsigned n) {
  struct ldata *v;
  struct cfg_item *w,*t;
  /* template sub-item array */
  for (v=p->item->first; v; v=v->l_next) {
    w = v->l_data;
    /* item copy with replaced first dependency mark */
    t = cfg_item_new();
    t->key = str_copy_new(w->key);
    str_subst_one_ui(&(t->key),d->dep_mark,id);
    t->optional = w->optional;
    t->dep_level = w->dep_level;
    t->data_type = w->data_type;
    t->data_rank = w->data_rank;
    t->data_ptr = cfg_data_pointer(w->data_ptr,w->data_type,
      w->data_rank+w->dep_level,id-1,n);
    t->data_dm1 = cfg_data_pointer(w->data_dm1,TYPE_UINT,
      w->data_rank+w->dep_level-1,id-1,n);
    t->data_dm2 = cfg_data_pointer(w->data_dm2,TYPE_UINT,
      w->data_rank+w->dep_level-2,id-1,n);
    t->read_fce = w->read_fce;
    /* recursive generation of sub-item tree */
    cfg_set_subsubitems(d,w,t,id,n);
    /* add sub-item to list */
    list_add_end_p(c->item,t);
    }
  }

/* define subitems of the given options according to saved template
  
   d - main config data struct
   p - parent item */
void cfg_set_subitems(struct cfg *d, struct cfg_item *p) {
  unsigned i,n,dl;
  struct ldata *s,*r;
  struct cfg_item *c,*t;
  struct queue *q;
  /* no sub-items in this options */
  if (!p->item->num)
    return;
  /* length of sub-item array */
  n = *((unsigned*)p->data_ptr);
  if (!n)
    return;
  q = queue_alloc();
  /* set of sub-items */
  for (s=p->item->first; s; s=s->l_next) {
    c = s->l_data;
    /* allocate data arrays in the whole subtree */
    queue_add(q,c);
    while (q->num) {
      t = queue_get(q);
      dl = str_char_count(t->key,d->dep_mark);
      cfg_alloc_data(t->data_ptr,t->data_type,t->data_rank+dl,n);
      cfg_alloc_data(t->data_dm1,TYPE_UINT,t->data_rank+dl-1,n);
      cfg_alloc_data(t->data_dm2,TYPE_UINT,t->data_rank+dl-2,n);
      for (r=t->item->first; r; r=r->l_next)
        queue_add(q,r->l_data);
      }
    /* sub-item array */
    dl = str_char_count(c->key,d->dep_mark);
    for (i=0; i<n; i++) {
      t = cfg_item_new();
      t->key = str_copy_new(c->key);
      str_subst_one_ui(&(t->key),d->dep_mark,i+1);
      t->optional = c->optional;
      t->data_type = c->data_type;
      t->data_rank = c->data_rank;
      t->data_ptr = cfg_data_pointer(c->data_ptr,c->data_type,
        c->data_rank+dl,i,n);
      t->data_dm1 = cfg_data_pointer(c->data_dm1,TYPE_UINT,
        c->data_rank+dl-1,i,n);
      t->data_dm2 = cfg_data_pointer(c->data_dm2,TYPE_UINT,
        c->data_rank+dl-2,i,n);
      t->read_fce = c->read_fce;
      /* sub-sub-items */
      cfg_set_subsubitems(d,c,t,i+1,n);
      /* add the sub-item to the main list */
      list_add_end_p(d->item,t);
      }
    }
  /* clean memory */
  queue_free(q);
  }

/* -------------------------------------------------------------------------- */
