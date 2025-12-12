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
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* add new data to the stack

   s    - pointer to the stack
   data - pointer to the data */
void stack_add(struct stack *s, void *data) {
   struct sdata *sp=NULL;
   sp = malloc(sizeof(struct sdata));
   if (!sp)
     msg_error("cannot allocate memory for stack cell",1);
   sp->s_data = data;
   if (!s->num)
     sp->s_next = NULL;
   else 
     sp->s_next = s->data;
   s->data = sp;
   (s->num)++;
   }

/* add integer to the stack and return actual number of items in the stack

   s   - pointer to the stack
   val - the integer */
long unsigned stack_iadd(struct stack *s, int val) {
  int *data=NULL;
  data = (int*)malloc(sizeof(int));
  if (!data)
    msg_error("cannot allocate memory for i-stack data",1);
  (*data) = val;
  stack_add(s,data);
  return(s->num);
  }

/* add unsigned integer to the stack and return actual number of items 

   s   - pointer to the stack
   val - the unsigned integer */
long unsigned stack_uadd(struct stack *s, unsigned val) {
  unsigned *data = NULL;
  data = (unsigned*)malloc(sizeof(unsigned));
  if (!data)
    msg_error("cannot allocate memory for u-stack data",1);
  (*data) = val;
  stack_add(s,data);
  return(s->num);
  }

/* add double real to the stack and return actual number of items in the stack

   s   - pointer to the stack
   val - the integer */
long unsigned stack_fadd(struct stack *s, double val) {
  double *data = NULL;
  data = (double*)malloc(sizeof(double));
  if (!data)
    msg_error("cannot allocate memory for f-stack data",1);
  (*data) = val;
  stack_add(s,data);
  return(s->num);
  }

/* add string to the stack and return actual number of items in the stack

   s   - pointer to the stack
   str - the string */
long unsigned stack_sadd(struct stack *s, char *str) {
  char *data = NULL;
  data = str_copy_new(str);
  stack_add(s,data);
  return(s->num);
  }

/* -------------------------------------------------------------------------- */
