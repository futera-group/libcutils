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
#include "cmn/arg.h"
#include "cmn/list.h"
#include "cmn/message.h"

/* -------------------------------------------------------------------------- */

/* print one command-line argument specification

   a     - command-line arguments
   d     - argument data structure
   space - specify if beginning space is printed or not */
void arg_print_item(struct arg *a, struct arg_dat *d, short space) {
  if (space)
    printf(" ");
  if (d->simple)
    printf("<%s>",d->desc);
  else if (d->optional) 
    printf("[-%c <%s>]",d->name[0],d->desc);
  else
    printf("-%c <%s>",d->name[0],d->desc);
  d->set=1;
  }

/* print list of bound command-line arguments

   a     - command-line arguments
   d     - argument data structure
   space - specify if beginning space is printed or not */
void arg_print_group(struct arg *a, struct arg_dat *d, short space) {
  unsigned n = 0;
  struct arg_dat *t;
  struct ldata *p;
  if (!d->bind) 
    arg_print_item(a,d,1);
  else {
    if (space)
      printf(" ");
    printf("{");
    for (p=a->data->first; p; p=p->l_next) {
      t = (struct arg_dat*)p->l_data;
      if (t->bind==d->bind)
        arg_print_item(a,t,n++);
      }
    printf("}");
    }
  }

/* print list of connected command-line arguments

   a  - command-line argument data
   id - first data structure in a group */
void arg_print_connected(struct arg *a, struct arg_dat *d) {
  short bracket;
  struct arg_dat *t,*x;
  struct ldata *p,*q;
  /* exit if the option is already printed */
  if (d->set)
    return;
  /* check if brackets are needed */
  bracket = 0;
  if (d->bind) {
    for (p=a->data->first; p; p=p->l_next) {
      t = (struct arg_dat*)p->l_data;
      if (t->bind==d->bind && t->excl) {
        bracket = 1;
        break;
        }
      }
    }
  /* print out items */
  if (bracket)
    printf(" {");
  /* group of connected items */
  arg_print_group(a,d,!bracket);
  if (d->bind)
    for (p=a->data->first; p; p=p->l_next) {
      t = (struct arg_dat*)p->l_data;
      if (t->bind==d->bind && t->excl)
        for (q=a->data->first; q; q=q->l_next) {
          x = (struct arg_dat*)q->l_data;
          if (t->excl==x && !x->set) {
            printf(" |");
            arg_print_connected(a,x);
            }
          }
      }
  /* group of excluded items */
  else if (d->excl && !d->excl->set) {
    printf(" |");
    arg_print_connected(a,d->excl);
    }
  if (bracket)
    printf("}");
  }

/* print usage of the program
   
   a       - list of possible and required arguments
   program - name of the program
   msg     - if not NULL, error message msg is printed below the usage
   quit    - exit program or not (0,1) */
void arg_print_usage(struct arg *a, char *program, char *msg, short quit) {
  short print_switch;
  struct arg_dat *d;
  struct ldata *p;
  /* initialization */
  if (a->data) {
    for (p=a->data->first; p; p=p->l_next)
      ((struct arg_dat*)p->l_data)->set = 0;
    }
  /* name of the program */
  printf("usage: %s",program);
  if (a->data) {
    /* switches */
    print_switch = 0;
    for (p=a->data->first; p; p=p->l_next)
      if (((struct arg_dat*)p->l_data)->swtch) {
        print_switch = 1;
        break;
        }
    if (print_switch) {
      printf(" [-");
      for (p=a->data->first; p; p=p->l_next) {
        d = (struct arg_dat*)p->l_data;
        if (d->swtch) {
          printf("%c",d->name[0]);
          d->set = 1;
          }
        }
      printf("]");
      }
    /* optional options */
    for (p=a->data->first; p; p=p->l_next) {
      d = (struct arg_dat*)p->l_data;
      if (!d->set && !d->simple && d->optional)
        arg_print_connected(a,d);
      }
    /* compulsory options */
    for (p=a->data->first; p; p=p->l_next) {
      d = (struct arg_dat*)p->l_data;
      if (!d->set && !d->simple)
        arg_print_connected(a,d);
      }
    /* simple arguments */
    for (p=a->data->first; p; p=p->l_next)
      if (!(((struct arg_dat*)p->l_data)->set))
        arg_print_connected(a,d);
    }
  printf("\n");
  /* error message */
  if (msg) 
    msg_error(msg,quit);
  }

/* -------------------------------------------------------------------------- */
