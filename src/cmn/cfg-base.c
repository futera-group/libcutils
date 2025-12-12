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
#include <stdlib.h>
#include "cmn/config.h"
#include "cmn/list.h"
#include "cmn/message.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for new config file struct */
struct cfg *cfg_new(void) {
  struct cfg *c = NULL;
  /* allocate memory */
  c = malloc(sizeof(struct cfg));
  if (!c)
    msg_error("cannot allocate memory for config file struct",1);
  /* initialization */
  c->verbose = 0;
  c->allow_free_lines = 1;
  c->allow_comment_lines = 1;
  c->allow_comment_in_lines = 1;
  c->ignore_unknown_key = 1;
  c->ignore_unknown_frmt = 0;
  c->comment_mark = '#';
  c->option_mark = '$';
  c->dep_mark = '*';
  c->cont_mark = '&';
  c->item = list_alloc();
  return(c);
  }

/* allocate memory for new config item data struct */
struct cfg_item *cfg_item_new(void) {
  struct cfg_item *c = NULL;
  /* allocate memory */
  c = malloc(sizeof(struct cfg_item));
  if (!c)
    msg_error("cannot allocate memory for config item data struct",1);
  /* initialization */
  c->key = NULL;
  c->dep_level = 0;
  c->optional = 0;
  c->assigned = 0;
  c->data_type = 0;
  c->data_rank = 0;
  c->data_ptr = NULL;
  c->data_dm1 = NULL;
  c->data_dm2 = NULL;
  c->item = list_alloc();
  c->read_fce = NULL;
  return(c);
  }

/* -------------------------------------------------------------------------- */

/* free memory of config-file item data saved in the list

   d - pointer to the item data */
void cfg_item_free(struct cfg_item *d) {
  if (d) {
    str_free(d->key);
    list_free(d->item,NULL);
    }
  }

/* free memory of config-file item data saved in the list - wrapper

   d - pointer to the item data */
void cfg_item_free_w(void *d) {
  cfg_item_free(d);
  }

/* free memory allocated for config file struct

   c - pointer to the config data */
void cfg_free(struct cfg *c) {
  if (c) {
    list_free(c->item,cfg_item_free_w);
    free(c);
    }
  }

/* -------------------------------------------------------------------------- */
