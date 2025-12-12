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
#include "cmn/config.h"
#include "cmn/message.h"
#include "cmn/string.h"
#include "cmn/types.h"

/* -------------------------------------------------------------------------- */

/* define one option in config file and connect it with a variable

   c    - pointer to config file setting struct
   key  - keyword expected in config file
   opt  - indicator if the keyword is optional or mandatory
   rank - rank of the data (scalar, array, ...)
   type - type of data connected with the keyword
   d    - pointer to data variable
   n1   - storage for rank1 data dimension
   n2   - storage for rank2 data dimension
   fce  - user-defined reading function */
void cfg_add_option(struct cfg *c, char *key, short opt, unsigned rank,
  short type, void *d, void *n1, void *n2,
  void(*fce)(void*,void*,void*,FILE*)) {
  struct cfg_item *p;
  if (!d)
    msg_error_f("attempt to set null data pointer to \"%s\" keyword",1,key);
  if (fce && type!=TYPE_UINT)
    msg_error_f("user-defined data keyword \"%s\" must be u-int type",1,key);
  p = cfg_item_new();
  p->key = str_copy_new(key);
  p->dep_level = 0;
  p->optional = opt;
  p->data_type = type;
  p->data_rank = rank;
  p->data_ptr = d;
  p->data_dm1 = n1;
  p->data_dm2 = n2;
  list_add_end_p(c->item,p);
  p->read_fce = fce;
  }

/* define dependent option in config file and connect it with a variable

   c    - pointer to config file setting struct
   par  - keyword of the parent option
   key  - keyword expected in config file
   opt  - indicator if the keyword is optional or mandatory
   rank - rank of the data (scalar, array, ...)
   type - type of data connected with the keyword
   d    - pointer to data variable
   n1   - storage for rank1 data dimension
   n2   - storage for rank2 data dimension,
   fce  - user-defined reading function */
void cfg_add_option_dep(struct cfg *c, char *par, char *key, short opt,
  unsigned rank, short type, void *d, void *n1, void *n2,
  void(*fce)(void*,void*,void*,FILE*)) {
  struct cfg_item *t,*p;
  if (!d)
    msg_error_f("attempt to set null data pointer to \"%s\" keyword",1,key);
  if (fce && type!=TYPE_UINT)
    msg_error_f("user-defined data keyword \"%s\" must be u-int type",1,key);
  /* parental item */
  t = cfg_check_keyword(c->item,par,1);
  if (!t)
    msg_error_f("parental keyword \"%s\" of \"%s\" is not defined",1,par,key);
  if (t->data_rank>0)
    msg_error_f("parental keyword \"%s\" of \"%s\" has non-zero data rank",1,
      par,key);
  if (t->data_type!=TYPE_UINT && 
      t->data_type!=TYPE_ULINT && 
      t->data_type!=TYPE_USINT)
    msg_error_f("parental keyword \"%s\" of \"%s\" does not have u-int"
      " data type",1,par,key);
  /* new item */
  p = cfg_item_new();
  p->key = str_copy_new(key);
  p->dep_level = str_char_count(p->key,c->dep_mark);
  if (p->dep_level!=(t->dep_level+1))
    msg_error_f("inconsistent dependency level of \"%s\" keyword and \"%s\""
      " parent",1,par,key);
  p->optional = opt;
  p->data_type = type;
  p->data_rank = rank;
  p->data_ptr = d;
  p->data_dm1 = n1;
  p->data_dm2 = n2;
  if (p->data_rank>=1 && !p->data_dm1)
    msg_error_f("storage for rank-1 data dimensions of keyword \"%s\" "
     "is not defined",1,key);
  if (!fce && p->data_rank>=2 && !p->data_dm2)
    msg_error_f("storage for rank-2 data dimensions of keyword \"%s\" "
     "is not defined",1,key);
  list_add_end_p(t->item,p);
  p->read_fce = fce;
  }

/* -------------------------------------------------------------------------- */
