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
#include "cmn/menu.h"
#include "cmn/message.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for menu data */
struct menu* menu_new(void) {
  struct menu *m = NULL;
  /* memory allocation */
  m = (struct menu*)malloc(sizeof(struct menu));
  if (!m)
    msg_error("cannot allocate memory for menu data",1);
  /* initialization */
  m->title = NULL;
  m->off_title = 0;
  m->off_item = 0;
  m->item = list_alloc();
  return(m);
  }

/* allocate and initialize menu data

   m     - pointer to menu struct
   ot    - offset for printing title 
   oi    - offset for printing items
   title - title of the menu */
struct menu* menu_new_init(unsigned ot, unsigned oi, char *title) {
  struct menu *m;
  /* memory allocation */
  m = menu_new();
  /* intialization */
  m->title = str_copy_new(title);
  m->off_title = ot;
  m->off_item = oi;
  return(m);
  }

/* allocate memory for menu item */
struct menu_item* menu_item_new(void) {
  struct menu_item *m = NULL;
  /* memory allocation */
  m = (struct menu_item*)malloc(sizeof(struct menu_item));
  if (!m)
    msg_error("cannot allocate memory for menu item",1);
  /* initialization */
  m->stop = 0;
  m->flag = NULL;
  m->desc = NULL;
  m->fce = NULL;
  m->fce_dat_1 = NULL;
  m->fce_dat_2 = NULL;
  return(m);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocate for menu item

   m - the menu item */
void menu_item_free(struct menu_item *m) {
  if (m) {
    str_free(m->flag);
    str_free(m->desc);
    free(m);
    }
  }

/* free memory allocate for menu item - wrapper

   m - the menu item */
void menu_item_free_w(void *m) {
  menu_item_free(m);
  }

/* free memory allocate for menu data

   m - the menu data */
void menu_free(struct menu* m) {
  if (m) {
    str_free(m->title);
    list_free(m->item,menu_item_free_w);
    }
  }

/* -------------------------------------------------------------------------- */
