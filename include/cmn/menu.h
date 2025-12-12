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

#ifndef ZF_LIB_CMN_MENU_H
#define ZF_LIB_CMN_MENU_H

#include "cmn/list.h"

/* -------------------------------------------------------------------------- */

/* menu item struct */
struct menu_item {
  short stop;                  /* specify if the item is terminating or not */
  char *flag;                  /* menu option */
  char *desc;                  /* menu item description */
  void (*fce)(void*,void*);    /* menu function */
  void *fce_dat_1;             /* pointer to menu function data */
  void *fce_dat_2;             /* pointer to menu function data */
  };

/* interactive menu struct */
struct menu {
  char *title;                          /* menu title */
  unsigned off_title;                   /* offset for printing title */
  unsigned off_item;                    /* offset for printing items */
  struct list *item;                    /* list of menu items */
  };

/* -------------------------------------------------------------------------- */

/* allocate memory for menu data */
struct menu* menu_new(void);
/* allocate and initialize menu data */
struct menu* menu_new_init(unsigned, unsigned, char*);
/* allocate memory for menu item */
struct menu_item* menu_item_new(void);

/* free memory allocate for menu data */
void menu_free(struct menu*);
/* free memory allocate for menu item */
void menu_item_free(struct menu_item*);

/* set one item of menu */
void menu_add_item(struct menu*, char*, char*, short,
  void (*fce)(void*,void*), void*, void*);

/* print out interactive menu */
void menu_print(struct menu*);

/* -------------------------------------------------------------------------- */

#endif
