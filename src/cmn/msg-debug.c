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

/* -------------------------------------------------------------------------- */

/* print debug message

   msg  - the message which will be printed 
   prnt - the message will be printed only if prnt!=0 */
void msg_debug(char *msg, short prnt) {
  if (prnt)
    printf("debug: %s\n",msg);
  }

/* print debug message and flusth the output stream
 
   msg - the message which will be printed out */
void msg_Debug(char *msg) {
  printf("DEBUG: %s\n",msg);
  fflush(stdout);
  }

/* pring function debug message and flush the output stream 

   fce - name of the function
   msg - the message which will be printed out */
void msg_Debug_fce(char *fce, char *msg) {
  printf("DEBUG: %s: %s\n",fce,msg);
  fflush(stdout);
  }

/* -------------------------------------------------------------------------- */
