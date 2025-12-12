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

#include "cmn/arg.h"
#include "cmn/list.h"
#include "cmn/message.h"

/* -------------------------------------------------------------------------- */

/* check compatibility of binding and exlusion setting

   a - command-line argument structure */
void arg_bind_check(struct arg *a) {
  struct arg_dat *t,*x;
  struct ldata *p,*q;
  if (a->data)
    for (p=a->data->first; p; p=p->l_next) {
      t = (struct arg_dat*)p->l_data;
      if (t->bind) 
        for (q=p->l_next; q; q=q->l_next) {
          x = (struct arg_dat*)q->l_data;
          if (t->bind==x->bind) {
            if (t->excl && t->excl->bind==t->bind)
              msg_error("binding and exclusion of the same options",1);
            if (t->optional!=x->optional) 
              msg_error("binding optional and mandatory option",1);
            }
          }
      }
  }

/* set option binding 

   a     - command-line argument structure
   d1,d2 - pointers to argument data to be bound */
void arg_bind(struct arg *a, struct arg_dat *d1, struct arg_dat *d2) {
  struct arg_dat *t;
  struct ldata *p;
  if (!d1->bind && !d2->bind) {
    d1->bind = d1;
    d2->bind = d1;
    }
  else if (!d1->bind && d2->bind)
    d1->bind = d2->bind;
  else if (d1->bind && !d2->bind)
    d2->bind = d1->bind;
  else if (d1->bind && d2->bind && d1->bind!=d2->bind) {
    for (p=a->data->first; p; p=p->l_next) {
      t = (struct arg_dat*)p->l_data;
      if (t->bind==d1->bind)
        t->bind = d2->bind;
      }
    }
  arg_bind_check(a);
  }

/* set option exclusion 

   a     - command-line argument structure
   d1,d2 - pointers to argument data to be bound */
void arg_exclude(struct arg *a, struct arg_dat *d1, struct arg_dat *d2) {
  d1->excl = d2;
  d2->excl = d1;
  arg_bind_check(a);
  }

/* -------------------------------------------------------------------------- */
