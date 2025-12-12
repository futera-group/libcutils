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
#include "cmn/message.h"

/* -------------------------------------------------------------------------- */

/* allocate new list structure */
struct list *list_alloc(void) {
  struct list *l = NULL;
  l = (struct list*)malloc(sizeof(struct list));
  if (!l)
    msg_error("cannnot allocate new list",1);
  l->num = 0;
  l->first = NULL;
  l->last = NULL;
  return(l);
  }

/* clean all the data saved in the list
 
   l - pointer to the list struct
   f - function for cleaning data */
void list_clean(struct list *l, void(f)(void*)) {
  struct ldata *lp,*lp2;
  if (!l) return;
  lp = l->first;
  while (lp) {
    lp2 = lp;
    lp = lp->l_next;
    if (lp2->l_data && f)
      f(lp2->l_data);
    free(lp2);
    }
  lp = l->first = NULL;
  l->num = 0;
  }

/* free memory allocated for the list
 
   l - pointer to the list struct
   f - function for cleaning data */
void list_free(struct list *l, void(f)(void*)) {
  list_clean(l,f);
  free(l);
  }

/* -------------------------------------------------------------------------- */
