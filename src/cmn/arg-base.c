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
#include "cmn/arg.h"
#include "cmn/list.h"
#include "cmn/message.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for the command-line argument structure */
struct arg* arg_new(void) {
  struct arg *a = NULL;
  /* memory allocation */
  a = (struct arg*)malloc(sizeof(struct arg));
  if (!a)
    msg_error("cannot allocate memory for command-line arguments",1);
  /* initialiation */
  a->help = 0;
  a->help_fce = NULL;
  a->data = list_alloc();
  return(a);
  }

/* allocate memory for argument data structures */
struct arg_dat* arg_dat_new(void) {
  struct arg_dat *a = NULL;
  /* memory allocation */
  a = (struct arg_dat*)malloc(sizeof(struct arg_dat));
  if (!a) 
    msg_error("cannot allocate memory for command-line argument data",1);
  /* initialization */
  a->name = NULL;
  a->desc = NULL;
  a->type = 0;
  a->swtch = 0;
  a->optional = 0;
  a->simple = 0;
  a->set = 0;
  a->bind = NULL;
  a->excl = NULL;
  a->val = NULL;
  return(a);
  }

/* -------------------------------------------------------------------------- */

/* wrapper of a node-cleaning function */
void arg_dat_free_w(void *d) {
  arg_dat_free(d);
  }

/* free memory allocated for a command-line argument structure */
struct arg* arg_free(struct arg *a) {
  if (a) {
    list_free(a->data,arg_dat_free_w);
    free(a);
    a = NULL;
    }
  return(a);
  }

/* free memory allocated for a command-line argument data structure */
struct arg_dat* arg_dat_free(struct arg_dat *d) {
  if (d) {
    str_free(d->name);
    str_free(d->desc);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
