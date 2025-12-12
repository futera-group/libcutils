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
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* print warning
  
   msg - warning message */
void msg_warn(char *msg) {
  printf("warning: %s\n",msg);
  }

/* print warning with user specified format
  
   msg - warning message */
void msg_warn_f(char *msg, ...) {
  char *t;
  va_list ap;
  va_start(ap,msg);
  t = str_merge_new("warning: ",msg);
  vprintf(t,ap);
  printf("\n");
  str_free(t);
  va_end(ap);
  }

/* -------------------------------------------------------------------------- */
