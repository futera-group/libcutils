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
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* check presence and compatibility
   
   a       - list of argument
   program - name of the program */
void arg_check(struct arg *a, char *program) {
  char *msg;
  struct arg_dat *t,*x;
  struct ldata *p,*q;
  if (a->data) {
    for (p=a->data->first; p; p=p->l_next) {
      t = (struct arg_dat*)p->l_data;
      /* check bound arguments */
      if (t->bind) {
        for (q=a->data->first; q; q=q->l_next) {
          x = (struct arg_dat*)q->l_data;
          if (t!=x && t->bind==x->bind && t->set!=x->set) {
            msg = str_new(45+str_length(t->desc)+str_length(x->desc));
            sprintf(msg,"arguments <%s> and <%s> must be set together",
              t->desc,x->desc);
            arg_print_usage(a,program,msg,1);
            }
          }
        }
      /* check excluded arguments */
      else if (t->excl) {
        if (t->set && t->excl->set) {
          msg = str_new(45+str_length(t->desc)+str_length(t->excl->desc));
          sprintf(msg,"arguments <%s> and <%s> must not be set together",
            t->desc,t->excl->desc);
          arg_print_usage(a,program,msg,1);
          }
        }
      /* check mandatory arguments */
      else if (!t->optional && !t->set) {
        if (!(t->excl && t->excl->set)) {
          msg = str_new(25+str_length(t->desc));
          sprintf(msg,"argument <%s> is not set",t->desc);
          arg_print_usage(a,program,msg,1);
          }
        }
      }
    }
  }

/* -------------------------------------------------------------------------- */
