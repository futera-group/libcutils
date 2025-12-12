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

#include "cmn/arg.h"
#include "cmn/list.h"
#include "cmn/message.h"
#include "cmn/string.h"
#include "cmn/types.h"

/* -------------------------------------------------------------------------- */

/* define switch and return its ID (offset in an array of arguments)
   
   a      - list of arguments
   name   - name of the switch
   single - 1 if switch is single operating (no other arguments are needed)
   val    - variable assigned to the switch */
struct arg_dat* arg_add_switch(struct arg *a, char *name,
  short single, short* val) {
  struct arg_dat *d;
  if (!name || name[0]=='\0')
    msg_error("name of argument switch is not defined",1);
  d = arg_dat_new();
  d->name = str_copy_new(name);
  d->type = TYPE_SINT;
  d->swtch = 1;
  d->optional = 1;
  d->simple = single;
  d->val = val;
  list_add_end_p(a->data,d);
  return(d);
  }

/* define option and return its ID (offset in an array of arguments)
   
   a    - list of arguments
   name - name of the option
   desc - short description of the option
   type - type of a variable assigned to the option
   opt  - optional (1) or compulsory (0) argument
   val  - the variable assigned to the option */
struct arg_dat* arg_add_option(struct arg *a, char *name, char *desc,
  short type, short opt, void *val) {
  struct arg_dat *d;
  if (!name || name[0]=='\0')
    msg_error("name of argument option is not defined",1);
  if (!desc || desc[0]=='\0')
    msg_error("description of argument option is not defined",1);
  d = arg_dat_new();
  d->name = str_copy_new(name);
  d->desc = str_copy_new(desc);
  d->type = type;
  d->optional = opt;
  d->val = val;
  list_add_end_p(a->data,d);
  return(d);
  }

/* define simple argument and return its ID (offset in an array of arguments)
   
   a    - list of arguments
   desc - short description of the argument
   type - type of a variable assigned to the argument
   opt  - optional (1) or compulsory (0) argument
   val  - the variable assigned to the argument */
struct arg_dat* arg_add_item(struct arg *a, char *desc, short type,
  short opt, void *val) {
  struct arg_dat *d;
  if (!desc || desc[0]=='\0')
    msg_error("name of argument item is not defined",1);
  d = arg_dat_new();
  d->desc = str_copy_new(desc);
  d->type = type;
  d->optional = opt;
  d->simple = 1;
  d->val = val;
  list_add_end_p(a->data,d);
  return(d);
  }

/* define simple argument and return its ID (offset in an array of arguments)
   
   a        - list of arguments
   fce_help - pointer to help function */
void arg_add_help(struct arg *a, void (*fce_help)(void)) {
  a->help = 1;
  a->help_fce = fce_help;
  }

/* -------------------------------------------------------------------------- */
