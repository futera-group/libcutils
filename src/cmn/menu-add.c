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

#include "cmn/list.h"
#include "cmn/menu.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* set one item of menu 

   m     - pointer to menu struct
   flag  - menu option
   desc  - menu line
   stop  - specify if the menu is printed again or not
   fce   - function assigned with the menu item
   d1,d2 - data for the fce function */
void menu_add_item(struct menu *m, char *flag, char *desc, short stop,
  void (*fce)(void*,void*), void *d1, void *d2) {
  struct menu_item *t;
  t = menu_item_new();
  t->stop = stop;
  t->flag = str_copy_new(flag);
  t->desc = str_copy_new(desc);
  t->fce = fce;
  t->fce_dat_1 = d1;
  t->fce_dat_2 = d2;
  list_add_end_p(m->item,t);
  }

/* -------------------------------------------------------------------------- */
