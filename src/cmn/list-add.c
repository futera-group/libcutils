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
#include <string.h>
#include "cmn/list.h"
#include "cmn/message.h"
#include "cmn/types.h"

/* -------------------------------------------------------------------------- */

/* add new data to begin of the list 
  
   l    - pointer to the list struct
   data - pointer to the data
   type - type of the data */
long unsigned list_add_begin(struct list *l, void *data, short type) {
  return(list_add_begin_p(l,type_copy(data,type)));
  }

/* add new data to begin of the list (pointer to data)
  
   l    - pointer to the list struct
   data - pointer to the data
   type - type of the data */
long unsigned list_add_begin_p(struct list *l, void *data) {
  struct ldata *lp = NULL;
  lp = malloc(sizeof(struct ldata));
  if (!lp)
    msg_error("cannot allocate memory for list cell",1);
  lp->l_data = data;
  lp->l_next = l->first;
  lp->l_prev = NULL;
  if (!l->num)
    l->first = l->last = lp;
  else {
    l->first->l_prev = lp;
    l->first = lp;
    }
  (l->num)++;
  return(l->num);
  }

/* add new data to begin of the list (user define type)
  
   l    - pointer to the list struct
   data - pointer to the data
   fcpy - function for copying user defined type data 
   size - size of the data */
long unsigned list_add_begin_t(struct list *l, void *data,
  void(fcpy)(void*,void*),
  unsigned size) {
  void *dt = NULL;
  dt = malloc(size);
  if (!dt)
    msg_error("cannot allocate memory for list data",1);
  fcpy(dt,data);
  return(list_add_begin_p(l,dt));
  }

/* -------------------------------------------------------------------------- */

/* add new data to end of the list 

   l    - pointer to the list struct
   data - pointer to the data
   type - type of the data */
long unsigned list_add_end(struct list *l, void *data, short type) {
  return(list_add_end_p(l,type_copy(data,type)));
  }

/* add new data to end of the list (pointer to data)

   l    - pointer to the list struct
   data - pointer to the data */
long unsigned list_add_end_p(struct list *l, void *data) {
  struct ldata *lp = NULL;
  lp = malloc(sizeof(struct ldata));
  if (!lp)
    msg_error("cannot allocate memory for list cell",1);
  lp->l_data = data;
  lp->l_next = NULL;
  lp->l_prev = l->last;
  if (!l->num)
    l->first = l->last=lp;
  else {
    l->last->l_next = lp;
    l->last = lp;
    }
  (l->num)++;
  return(l->num);
  }

/* add new data to end of the list (user define type)
  
   l    - pointer to the list struct
   data - pointer to the data
   fcpy - function for copying user defined type data 
   size - data size in bytes */
long unsigned list_add_end_t(struct list *l, void *data,
  void(fcpy)(void*,void*),
  unsigned size) {
  void *dt = NULL;
  dt = malloc(size);
  if (!dt)
    msg_error("cannot allocate memory for list data",1);
  fcpy(dt,data);
  return(list_add_end_p(l,dt));
  }

/* -------------------------------------------------------------------------- */

/* add new data to specific position in the list 

   l    - pointer to the list struct
   data - pointer to the data
   type - type of the data 
   id   - position in the list */
long unsigned list_add_pos(struct list *l, void *data, short type,
  long unsigned id) {
  return(list_add_pos_p(l,type_copy(data,type),id));
  }

/* add new data to specific position in the list (pointer to data)

   l    - pointer to the list struct
   data - pointer to the data
   id   - position in the list */
long unsigned list_add_pos_p(struct list *l, void *data, long unsigned id) {
  unsigned i;
  struct ldata *ln,*lp;
  ln = NULL;
  lp = l->first;
  if (id==0)
    list_add_begin_p(l,data);
  else if (id==l->num)
    list_add_end_p(l,data);
  else if (id<l->num) {
    for (i=0; i<id; i++)
      lp = lp->l_next;
    ln = malloc(sizeof(struct ldata));
    if (!ln)
      msg_error("cannot allocate memory for list cell",1);
    lp->l_data = data;
    ln->l_next = lp;
    ln->l_prev = lp->l_prev;
    ln->l_prev->l_next = ln;
    ln->l_next->l_prev = ln;
    (l->num)++;
    }
  return(l->num);
  }

/* add new data to specific position in the list (user defined type)
  
   l    - pointer to the list struct
   data - pointer to the data
   fcpy - function for copying user defined type data
   size - size of the data
   id   - specific position in the list */
long unsigned list_add_pos_t(struct list *l, void *data,
  void(fcpy)(void*,void*), unsigned size, long unsigned id) {
  void *dt = NULL;
  dt = malloc(size);
  if (!dt)
    msg_error("cannot allocate memory for list data",1);
  fcpy(dt,data);
  return(list_add_pos_p(l,dt,id));
  }

/* -------------------------------------------------------------------------- */
