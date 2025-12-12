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

#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "cmn/file.h"
#include "cmn/message.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* detect if given directory exists

   d - path to the directory */
int dir_exist(char *d) {
  struct stat st;
  if (stat(d,&st)==0 && S_ISDIR(st.st_mode))
   return(1);
  return(0);
  }

/* get the current directory path
 
   s - storage for the path */
void dir_current(char *s) {
  if (!getcwd(s,str_length(s)))
    msg_error("cannot get name of current directory",1);
  }

/* return new string with current working directory path
 
   s - storage for the path */
char* dir_current_new(void) {
  char *s;
  s = getcwd(NULL,MAXPATHLEN);
  if (!s)
    msg_error("cannot get name of current directory",1);
  return(s);
  }

/* -------------------------------------------------------------------------- */

/* create directory with given mask
 
   path - directory path
   mode - mask
   code - exit if not 0 when an error occurs */
int dir_make_one(char *path, mode_t mode, short code) {
  char msg[80];
  if (mkdir(path,mode)<0) {
    switch (errno) {
      case EACCES:       sprintf(msg,"permission denied");            break;
      case EDQUOT:       sprintf(msg,"disk quota exceeded");          break;
      case EEXIST:       sprintf(msg,"directory already exists");     break;
      case EFAULT:       sprintf(msg,"bad address space");            break;
      case EIO:          sprintf(msg,"I/O error");                    break;
      case ELOOP:        sprintf(msg,"too many symbolic links");      break;
      case EMLINK:       sprintf(msg,"max number of links exceeded"); break;
      case ENAMETOOLONG: sprintf(msg,"too long path");                break;
      case ENOENT:       sprintf(msg,"no valid entry");               break;
      case ENOSPC:       sprintf(msg,"no space on the file system");  break;
      case ENOTDIR:      sprintf(msg,"not a directory");              break;
      case EROFS:        sprintf(msg,"read-only file system");        break;
      default:           sprintf(msg,"unknown error");                break;
      }
    msg_error_f("cannot make directory %s (%s)",code,path,msg);
    return(0);
    }
  return(1);
  }

/* create directory with given mask
 
   path - directory path
   mode - mask
   code - exit if not 0 when an error occurs */
int dir_make(char *path, mode_t mode, short code) {
  char **s,*d0,*d1;
  unsigned i,n;
  int ok = 1;
  /* split to sub-directories */
  s = str_split(path,'/',&n);
  for (i=0,d0=NULL,d1=NULL; i<n; i++) {
    /* root directory */
    if (i==0) {
      if (path[0]=='/')
        d1 = str_merge_new("/",s[i]);
      else
        d1 = str_copy_new(s[i]);
      }
    /* sub-directories */
    else
      d1 = str_file_path_new(d0,s[i]);
    /* check path */
    if (!file_exist(d1)) 
      ok = dir_make_one(d1,mode,code);
    str_free(d0);
    d0 = d1;
    d1 = NULL;
    if (!ok)
      break;
    }
  str_free(d0);
  return(ok);
  }

/* -------------------------------------------------------------------------- */

/* change current working directory

   path - path to the directory
   code - exit if not 0 when an error occurs */
int dir_change(char *path, short code) {
  char msg[80];
  if (chdir(path)<0) {
    switch (errno) {
      case EFAULT:       sprintf(msg,"bad address space");          break;
      case ENAMETOOLONG: sprintf(msg,"too long path");              break;
      case ENOENT:       sprintf(msg,"no such directory");          break;
      case ENOMEM:       sprintf(msg,"insufficient kernel memory"); break;
      case ENOTDIR:      sprintf(msg,"not a directory");            break;
      case EACCES:       sprintf(msg,"permission denied");          break;
      case ELOOP:        sprintf(msg,"too many symbolic links");    break;
      case EIO:          sprintf(msg,"I/O error");                  break;
      default:           sprintf(msg,"unknown error");              break;
      }
    msg_error_f("cannot change to directory %s (%s)",1,path,msg);
    return(0);
    }
  return(1);
  }

/* -------------------------------------------------------------------------- */

/* open directory and return its pointer

   path - path to the directory
   code - exit if not 0 when an error occurs */
DIR *dir_open(char *path, short code) {
  char msg[80];
  DIR *dir;
  dir = opendir(path);
  if (!dir) {
    switch (errno) {
      case EACCES:  sprintf(msg,"permission denied");              break;
      case EMFILE:  sprintf(msg,"too many used file descriptors"); break;
      case ENFILE:  sprintf(msg,"too many opened file");           break;
      case ENOENT:  sprintf(msg,"no such directory");              break;
      case ENOMEM:  sprintf(msg,"insufficient memory");            break;
      case ENOTDIR: sprintf(msg,"not a directory");                break;
      default:      sprintf(msg,"unknown error");                  break;
      }
    msg_error_f("cannot open directory %s (%s)",1,path,msg);
    }
  return(dir);
  }

/* open the current working directory for reading */
DIR *dir_open_current(void) {
  char *dir_name;
  DIR *dir;
  dir_name = dir_current_new();
  dir = dir_open(dir_name,1);
  str_free(dir_name);
  return(dir);
  }

/* close opened directory

   dir - pointer to the directory */
void dir_close(DIR *dir) {
  closedir(dir);
  }

/* -------------------------------------------------------------------------- */
