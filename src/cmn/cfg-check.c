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
#include <string.h>
#include "cmn/config.h"
#include "cmn/message.h"

/* -------------------------------------------------------------------------- */

/* check whether all mandatory options are assigned

   c - pointer to config file setting */
void cfg_check_assignment(struct cfg *c) {
  struct ldata *p;
  struct cfg_item *d;
  for (p=c->item->first; p; p=p->l_next) {
    d = p->l_data;
    if (!d->optional && !d->assigned)
      msg_error_f("keyword \"%s\" was not found",1,d->key);
    }
  }

/* search the item list and find element with specified keyword

   t   - the item list
   key - the keyword
   rec - recursive search */
struct cfg_item* cfg_check_keyword(struct list *t, char *key, short rec) {
  struct cfg_item *c,*r = NULL;
  struct ldata *p;
  /* search the item list */
  for (p=t->first; p; p=p->l_next) {
    c = p->l_data;
    if (strcmp(c->key,key)==0) {
      r = c;
      break;
      }
    /* recursion */
    if (rec) {
      r = cfg_check_keyword(c->item,key,rec);
      if (r) 
        break;
      }
    }
  return(r);
  }

/* -------------------------------------------------------------------------- */
