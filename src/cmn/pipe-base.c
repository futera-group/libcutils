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
#include <stdio.h>
#include "cmn/message.h"

/* -------------------------------------------------------------------------- */

/* open pipe for given command and return its stream pointer

   cmd  - the command line
   mode - opening mode (e.i. "r","w","rw",...) */
FILE *pipe_open(char *cmd, char *mode) {
  FILE *f = NULL;
  f = popen(cmd,mode);
  if (!f)
    msg_error_f("cannot open pipe for \"%s\"",1,cmd);
  return(f);
  }

/* close open pipe stream
 
   f - the stream */
void pipe_close(FILE *f) {
  pclose(f);
  }

/* -------------------------------------------------------------------------- */
