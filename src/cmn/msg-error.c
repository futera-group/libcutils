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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* error function struct */
struct {
  short ef_f;
  void (*ef_fce)(void*, void*, void*);
  void *ef_a1;
  void *ef_a2;
  void *ef_a3;
  } msg_ef = {1,NULL,NULL,NULL,NULL};

/* set pointer to error function

   e_fce - pointer to error function
   e_a1  - first argument to error function 
   e_e2  - second argument to the error function */
void msg_error_fce_set(void(*e_fce)(void*,void*,void*),
  void *e_a1, void *e_a2, void *e_a3) {
  msg_ef.ef_fce = e_fce;
  msg_ef.ef_a1 = e_a1;
  msg_ef.ef_a2 = e_a2;
  msg_ef.ef_a3 = e_a3;
  }

/* print error message
   
   msg  - error message
   quit - exit program if 1 */
void msg_error(char *msg, short quit) {
  /* print error message */
  printf("error: %s\n",msg);
  /* pre-defined error function */
  if (msg_ef.ef_fce!=NULL && msg_ef.ef_f) {
    msg_ef.ef_f = 0;
    msg_ef.ef_fce(msg_ef.ef_a1,msg_ef.ef_a2,msg_ef.ef_a3);
    }
  /* quit */
  if (quit)
    exit(1);
  msg_ef.ef_f = 1;
  }

/* print formatted error message
   
   msg  - error message
   quit - exit program if 1 */
void msg_error_f(char *msg, short quit, ...) {
  static char *t;
  va_list ap;
  /* print error message */
  va_start(ap,quit);
  t = str_merge_new("error: ",msg);
  vprintf(t,ap);
  printf("\n");
  str_free(t);
  va_end(ap);
  /* pre-defined error function */
  if (msg_ef.ef_fce!=NULL && msg_ef.ef_f) {
    msg_ef.ef_f = 0;
    msg_ef.ef_fce(msg_ef.ef_a1,msg_ef.ef_a2,msg_ef.ef_a3);
    }
  /* quit */
  if (quit) exit(1);
  msg_ef.ef_f = 1;
  }

/* -------------------------------------------------------------------------- */
