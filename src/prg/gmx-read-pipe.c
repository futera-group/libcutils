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
#include <cmn/file.h>
#include <cmn/pipe.h>
#include <cmn/string.h>

/* -------------------------------------------------------------------------- */

/* Call pre-processor on given file and send data to open pipe
 
   name   - name of the file
   def    - array of pre-processor definitions
   n_defs - number of preprocessor definitions */
FILE *gmx_read_pipe_open(char *name, char **def, unsigned n_defs) {
  char *prep,*cmd,*str;
  unsigned i;
  FILE *f;
  /* pre-processing command */
  prep = file_exec_path_new("cpp");
  if (n_defs) {
    cmd = str_sprintf("%s -P -w ",prep);
    for (i=0; i<n_defs; i++) {
      str = str_sprintf("-D%s ",def[i]);
      str_append(&cmd,str);
      }
    str_append(&cmd,name);
    }
  else
    cmd = str_sprintf("%s -P -w %s",prep,name);
  prep = str_free(prep);
  /* send output to pipe */
  f = pipe_open(cmd,"r");
  return(f);
  }

/* -------------------------------------------------------------------------- */
