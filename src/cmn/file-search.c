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

/* -------------------------------------------------------------------------- */

#include <string.h>
#include <sys/stat.h>
#include "cmn/dir.h"
#include "cmn/file.h"
#include "cmn/list.h"
#include "cmn/queue.h"
#include "cmn/string.h"

/* search and create list of all file of specified type

   root_dir - name of directory to search
   rec   - recursive search through all subdirectories
   ext   - file extension */
struct list *file_search(char *root_dir, short rec, char *ext) {
  char *dir_name,*file_name,*file_ext;
  DIR *dir;
  struct dirent *dir_info;
  struct stat f_stat;
  struct queue *q_dir;
  struct list *t;
  t = list_alloc();
  q_dir = queue_alloc();
  /* initial directory */
  if (root_dir && root_dir[0])
    queue_sadd(q_dir,root_dir);
  else
    queue_sadd(q_dir,".");
  /* search loop */
  while (q_dir->num) {
    dir_name = queue_get(q_dir);
    /* open directory */
    if (dir_name && dir_exist(dir_name)) {
      dir = dir_open(dir_name,1);
      dir_info = readdir(dir);
      /* list all files */
      while (dir_info) {
        file_name = str_file_path_new(dir_name,dir_info->d_name);
        if (file_name) {
          file_stat(file_name,&f_stat,1);
          /* regular file or link */
          if (S_ISREG(f_stat.st_mode) || S_ISLNK(f_stat.st_mode)) {
            file_ext = str_file_ext_new(dir_info->d_name);
            if (file_ext && strcmp(file_ext,ext)==0)
              list_add_end_p(t,file_dat_new(file_name));
            str_free(file_ext);
            }
          /* recursive sub-directory search */
          else if (rec && S_ISDIR(f_stat.st_mode) &&
            strcmp(dir_info->d_name,".")!=0 &&
            strcmp(dir_info->d_name,"..")!=0) {
            queue_sadd(q_dir,file_name);
            }
          dir_info = readdir(dir);
          str_free(file_name);
          }
        }
      /* close directory */
      dir_close(dir);
      }
    str_free(dir_name);
    }
  /* clean memory */
  queue_free(q_dir);
  return(t);
  }

/* -------------------------------------------------------------------------- */
