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
#include <string.h>
#include "cmn/stack.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* return pointer to the first data item in the stack

   s - poitner to the stack */
void *stack_get(struct stack *s) {
  void *data = NULL;
  struct sdata *sp;
  if (s->num>0) {
    sp = s->data;
    s->data = sp->s_next;
    data = sp->s_data;
    free(sp);
    (s->num)--;
    }
  return(data);
  }

/* get integer from the stack and return actual number of items in the stack

   s   - pointer to the stack
   val - the integer */
long unsigned stack_iget(struct stack *s, int *val) {
  void *data;
  data = stack_get(s);
  (*val) = (*(int*)data);
  free(data);
  return(s->num);
  }

/* get unsigned integer from the stack and return actual number of items

   s   - pointer to the stack
   val - the unsigned integer */
long unsigned stack_uget(struct stack *s, unsigned *val) {
  void *data;
  data = stack_get(s);
  (*val) = (*(unsigned*)data);
  free(data);
  return(s->num);
  }

/* get double real from the stack and return actual number of items

   s   - pointer to the stack
   val - the unsigned integer */
long unsigned stack_fget(struct stack *s, double *val) {
  void *data;
  data = stack_get(s);
  (*val) = (*(double*)data);
  free(data);
  return(s->num);
  }

/* get string from the stack and return actual number of items in the stack

   s   - pointer to the stack
   str - the string */
long unsigned stack_sget(struct stack *s, char *str) {
  void *data;
  data = stack_get(s);
  strcpy(str,(char*)data);
  free(data);
  return(s->num);
  }

/* -------------------------------------------------------------------------- */
