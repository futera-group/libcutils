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

#ifndef ZF_LIB_CMN_LIST_H
#define ZF_LIB_CMN_LIST_H

/* -------------------------------------------------------------------------- */

/* structs and types */
struct ldata {
  struct ldata *l_next;
  struct ldata *l_prev;
  void *l_data;
  };

struct list {
  long unsigned num;
  struct ldata *first;
  struct ldata *last;
  };

/* -------------------------------------------------------------------------- */

/* functions for allocating and destroying lists */

/* allocate new list structure */
struct list *list_alloc(void);
/* clean all the data saved in the list */
void list_clean(struct list*, void(f)(void*));
/* free memory allocated for the list */
void list_free(struct list*, void(f)(void*));

/* function for adding new data */

/* add new data to begin of the list */
long unsigned list_add_begin(struct list*, void*, short);
/* add new data to begin of the list (pointer to data) */
long unsigned list_add_begin_p(struct list*, void*);
/* add new data to begin of the list (user defined type) */
long unsigned list_add_begin_t(struct list*, void*, void(fcpy)(void*,void*),
  unsigned);
/* add new data to end of the list */
long unsigned list_add_end(struct list*, void*, short);
/* add new data to end of the list (pointer to data) */
long unsigned list_add_end_p(struct list*, void*);
/* add new data to end of the list (user defined type) */
long unsigned list_add_end_t(struct list*, void*, void(fcpy)(void*,void*),
  unsigned);
/* add new data to specific position in the list */
long unsigned list_add_pos(struct list*, void*, short, long unsigned);
/* add new data to specific position in the list (pointer to data) */
long unsigned list_add_pos_p(struct list*, void*, long unsigned);
/* add new data to specific position in the list (user defined type) */
long unsigned list_add_pos_t(struct list*, void*, void(fcpy)(void*,void*),
  unsigned, long unsigned);

/* functions for deleting data */

/* delete data from the specific pointer in the list */
long unsigned list_del_pointer(struct list*, struct ldata*, void(f)(void*));
/* delete data from the begin of the list */
long unsigned list_del_begin(struct list*, void(f)(void*));
/* delete data from the end of the list */
long unsigned list_del_end(struct list*, void(f)(void*));
/* delete data from the specific position in the list */
long unsigned list_del_pos(struct list*, unsigned, void(f)(void*));

/* function for getting data */

/* get data from the begin of the list */
void *list_get_begin(struct list*);
/* get data from the end of the list */
void *list_get_end(struct list*);
/* get data from the specific position of the list */
void *list_get_pos(struct list*, long unsigned);

/* functions for searching data list nodes  */

/* find element in list according to specified comparing function */
long unsigned list_find_id_t(struct list*, void*, short(f)(void*,void*));
/* find element in list according to specified comparing function */
void *list_find_data_t(struct list*, void*, short(f)(void*,void*));

/* functions for sorting list nodes */

/* reorder elements in the list according to specified element IDs */
void list_sort_id(struct list*, unsigned*);
/* sort the elements in the list by provided comparing function */
void list_sort(struct list*, short(f)(void*,void*));

/* conversion to different containers */

/* convert data saved in a list to vector */
void *list_conv_vec(struct list*, short, long unsigned*);

/* -------------------------------------------------------------------------- */

#endif
