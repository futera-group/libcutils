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
#include "cmn/list.h"
#include "cmn/menu.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* read menu item selection and perform appropriate function

   m - pointer to menu struct */
short menu_print_ask(struct menu *m) {
  char *answer;
  short repeat,stop;
  struct ldata *p;
  struct menu_item *t;
  /* initialization */
  repeat = 1;
  stop = 0;
  /* read answer */
  printf("\n");
  while (repeat) {
    if (m->off_item)
      printf("%*cAnswer: ",m->off_item,' ');
    else
      printf("Answer: ");
    answer = str_read_word_new();
    /* compare with menu items */
    for (p=m->item->first; p; p=p->l_next) {
      t = (struct menu_item*)p->l_data;
      if (str_compare(t->flag,answer)) {
        /* call menu function */
        if (t->fce) {
          printf("\n");
          t->fce(t->fce_dat_1,t->fce_dat_2);
          }
        stop = t->stop;
        repeat = 0;
        break;
        }
      }
    }
  return(stop);
  }

/* print out interactive menu

   m - pointer to menu struct */
void menu_print(struct menu *m) {
  unsigned n,flag_length;
  struct ldata *p;
  struct menu_item *t;
  /* set length of menu flag */
  flag_length = 0;
  for (p=m->item->first; p; p=p->l_next) {
    t = (struct menu_item*)p->l_data;
    n = str_length(t->flag);
    if (n>flag_length)
      flag_length = n;
    }
  /* print menu */
  for (;;) {
    /* title */
    if (m->off_title)
      printf("\n%*c%s\n\n",m->off_title,' ',m->title);
    else
      printf("\n%s\n\n",m->title);
    /* items */
    for (p=m->item->first; p; p=p->l_next) {
      t = (struct menu_item*)p->l_data;
      if (m->off_item)
        printf("%*c%*s: %s\n",m->off_item,' ',flag_length,t->flag,t->desc);
      else
        printf("%*s: %s\n",flag_length,t->flag,t->desc);
      }
    if (menu_print_ask(m))
      break;
    }
  }

/* -------------------------------------------------------------------------- */
