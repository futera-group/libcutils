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
#include "cmn/list.h"

/* -------------------------------------------------------------------------- */

/* delete data from the specific pointer in the list
 
   l - pointer to the list struct 
   p - the pointer to one record in the list
   f - memory cleaning function */
long unsigned list_del_pointer(struct list *l, struct ldata *p,
  void(f)(void*)) {
  /* shift pointers to next and previous record */
  if (p==l->first) {
    l->first = l->first->l_next;
    if (l->first)
      l->first->l_prev = NULL;
    }
  else if (p==l->last) {
    l->last = l->last->l_prev;
    if (l->last)
      l->last->l_next = NULL;
    }
  else {
    p->l_prev->l_next = p->l_next;
    p->l_next->l_prev = p->l_prev;
    }
  /* delete specified record */
  if (p->l_data && f) {
    f(p->l_data);
    p->l_data = NULL;
    }
  free(p);
  /* decrease number of records */
  (l->num)--;
  return(l->num);
  }

/* delete data from the begin of the list 
 
   l - pointer to the list struct
   f - memory cleaning function */
long unsigned list_del_begin(struct list *l, void(f)(void*)) {
  return(list_del_pointer(l,l->first,f));
  }

/* delete data from the end of the list 
 
   l - pointer to the list struct
   f - memory cleaning function */
long unsigned list_del_end(struct list *l, void(f)(void*)) {
  return(list_del_pointer(l,l->last,f));
  }

/* delete data from the specific position in the list 
 
   l  - pointer to the list struct 
   id - the position in the list
   f  - memory cleaning function */
long unsigned list_del_pos(struct list *l, unsigned id, void(f)(void*)) {
  unsigned i;
  struct ldata *lp = l->first;
  if (id==0)
    return(list_del_pointer(l,l->first,f));
  else if (id==l->num)
    return(list_del_pointer(l,l->last,f));
  else if (id<l->num) {
    for (i=0; i<id; i++)
      lp=lp->l_next;
    return(list_del_pointer(l,lp,f));
    }
  return(l->num);
  }

/* -------------------------------------------------------------------------- */
