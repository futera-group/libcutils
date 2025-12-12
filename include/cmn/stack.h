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

#ifndef ZF_LIB_CMN_STACK_H
#define ZF_LIB_CMN_STACK_H

/* -------------------------------------------------------------------------- */

/* struct and types */
struct sdata {
  struct sdata *s_next;
  void *s_data;
  };

struct stack {
  long unsigned num;
  struct sdata *data;
  };

/* -------------------------------------------------------------------------- */

/* allocate new stack structure */
struct stack *stack_alloc(void);
/* destroy allocated stack structure */
void stack_free(struct stack*);

/* add new unspecified dat to the stack */
void stack_add(struct stack*, void*);
/* add integer to the stack */
long unsigned stack_iadd(struct stack*, int);
/* add unsigned integer to the stack */
long unsigned stack_uadd(struct stack*, unsigned);
/* add double real to the stack */
long unsigned stack_fadd(struct stack*, double);
/* add string to the stack */
long unsigned stack_sadd(struct stack*, char*);

/* return pointer to the highest item in the stack */
void *stack_get(struct stack*);
/* get integer from the stack */
long unsigned stack_iget(struct stack*, int*);
/* get unsigned integer from the stack */
long unsigned stack_uget(struct stack*, unsigned*);
/* get double real from the stack */
long unsigned stack_fget(struct stack*, double*);
/* get string from the stack */
long unsigned stack_sget(struct stack*, char*);

/* -------------------------------------------------------------------------- */

#endif
