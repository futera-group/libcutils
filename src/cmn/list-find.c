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

/* find element in list according to specified comparing function 

   l - pointer to list struct
   d - data for comparing function
   f - comparing function */
long unsigned list_find_id_t(struct list *l, void *d, short(f)(void*,void*)) {
  unsigned n = 0;
  struct ldata *lp = l->first;
  if (!l) return(l->num);
  while (lp) {
    if (f(lp->l_data,d))
      return(n);
    lp = lp->l_next;
    n++;
    }
  return(l->num);
  }

/* find element in list according to specified comparing function

   l - pointer to list struct
   d - data for comparing function
   f - comparing function */
void *list_find_data_t(struct list *l, void *d, short(f)(void*,void*)) {
  struct ldata *lp = l->first;
  if (!l) return(NULL);
  while (lp) {
    if (f(lp->l_data,d))
      return(lp->l_data);
    lp = lp->l_next;
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
