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
#include "cmn/message.h"
#include "cmn/stack.h"

/* -------------------------------------------------------------------------- */

/* allocate new stack structure */
struct stack *stack_alloc(void) {
  struct stack *s;
  s = (struct stack*)malloc(sizeof(struct stack));
  if (!s)
    msg_error("cannot allocate new stack",1);
  s->data = NULL;
  s->num = 0;
  return(s);
  }

/* destroy the allocated stack

   s - pointer to the stack */
void stack_free(struct stack *s) {
  struct sdata *sp,*sp2;
  sp = s->data;
  while (sp) {
    sp2 = sp;
    sp = sp->s_next;
    free(sp2->s_data);
    free(sp2);
    }
  free(s);
  }

/* -------------------------------------------------------------------------- */
