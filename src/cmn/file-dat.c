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
#include "cmn/file.h"
#include "cmn/message.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* create new file info data struct from given path 
  
   p - patch to the file */
struct file_dat *file_dat_new(char *p) {
  struct file_dat *f = NULL;
  /* allocate memory */
  f = (struct file_dat*)malloc(sizeof(struct file_dat));
  if (!f)
    msg_error("cannot allocate memory for new file info data struct",1);
  /* data assignment */
  if (p) {
    /* directory */
    f->path = str_file_dir_new(p);
    if (!f->path || !f->path[0])
      f->path = str_copy_new(".");
    /* file name */
    f->name = str_file_name_new(p);
    }
  else {
    f->path = NULL;
    f->name = NULL;
    }
  return(f);
  }

/* free memory allocated for the file data struct
 
   f - the file data struct */
void file_dat_free(struct file_dat *f) {
  if (f) {
    str_free(f->name);
    str_free(f->path);
    free(f);
    }
  }

/* -------------------------------------------------------------------------- */
