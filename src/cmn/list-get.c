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

/* get data from the begin of the list 
 
   l - pointer to the list struct */
void *list_get_begin(struct list *l) {
  void *d = NULL;
  if (l->num)
    d=l->first->l_data;
  return(d);
  }

/* get data from the end of the list 
 
   l - pointer to the list struct */
void *list_get_end(struct list *l) {
  void *d = NULL;
  if (l->num)
    d=l->last->l_data;
  return(d);
  }

/* get data from the specific position of the list 
 
   l  - pointer to the list struct 
   id - the position in the list */
void *list_get_pos(struct list *l, long unsigned id) {
  unsigned i;
  struct ldata *lp = l->first;
  void *d = NULL;
  if (id==0)
    d = list_get_begin(l);
  else if (id==l->num)
    d = list_get_end(l);
  else if (id<l->num) {
    for (i=0; i<id; i++)
      lp = lp->l_next;
    d = lp->l_data;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
