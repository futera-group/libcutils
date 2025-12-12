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
#include "cmn/list.h"
#include "cmn/message.h"
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* reorder elements in the list according to specified element IDs

   t - pointer to the list struct
   r - reordering array with element IDs */
void list_sort_id(struct list *t, unsigned *r) {
  short found;
  unsigned i,j;
  long unsigned n;
  struct ldata *p,**s;
  /* save data pointers */
  n = t->num;
  s = vec_talloc(sizeof(struct ldata*),n);
  for (p=t->first,i=0; p; p=p->l_next,i++)
    s[i] = p->l_data;
  /* recreate the array in given order */
  list_clean(t,NULL);
  for (i=0; i<n; i++) {
    for (j=0,found=0; j<n; j++)
      if (r[i]==j) {
        list_add_end_p(t,s[j]);
        found = 1;
        break;
        }
    if (!found)
      msg_error("missing id in list-reordering array",1);
    }
  /* clean memory */
  vec_tfree(s);
  }

/* -------------------------------------------------------------------------- */

/* sort the elements in the list by provided comparing function

   s - pointer to list struct
   f - comparing function */
void list_sort(struct list *s, short(f)(void*,void*)) {
  short change;
  void *d;
  struct ldata *s1,*s2;
  if (s->num<=1)
    return;
  /* bubble sort - repeating loop */
  for (s1=s->first; s1; s1=s1->l_next) {
    change=0;
    /* bubble sort - sorting loop */
    s2 = s->last;
    for (s2=s->last; s2!=s1; s2=s2->l_prev) {
      if (f((s2->l_prev)->l_data,s2->l_data)>0) {
        d = (s2->l_prev)->l_data;
        (s2->l_prev)->l_data = s2->l_data;
        s2->l_data = d;
        change++;
        }
      }
    /* check state */
    if (!change)
      break;
    }
  }

/* -------------------------------------------------------------------------- */
