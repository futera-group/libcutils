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
#include <string.h>
#include <unistd.h>
#include "cmn/dir.h"
#include "cmn/message.h"
#include "cmn/string.h"
#include "cmn/vector.h"

#define PATH_ABSOLUTE 1
#define PATH_RELATIVE 2
#define PATH_NAMEONLY 3

/* -------------------------------------------------------------------------- */

/* check path to file and return its type

   path - the given path */
short file_exec_path_type(char *path) {
  if (path[0]=='/')
    return(PATH_ABSOLUTE);
  else if (path[0]=='.' && path[1]=='/')
    return(PATH_RELATIVE);
  else if (path[0]=='.' && path[1]=='.' && path[2]=='/')
    return(PATH_RELATIVE);
  return(PATH_NAMEONLY);
  }

/* check if the given file exist and is executable

   path - path to the given file */
void file_exec_path_check(char *path) {
  if (access(path,F_OK)!=0)
    msg_error_f("file %s doesn't exist",1,path);
  if (access(path,X_OK)!=0)
    msg_error_f("file %s is not executable",1,path);
  }

/* find path to the given executable program

   prog - name of the program */
void file_exec_path_find(char *prog) {
  short found;
  unsigned i,n;
  char *path,**pvec,*p;
  /* retrieve path variable */
  path = getenv("PATH");
  if (!path)
    msg_error("environment variable PATH is not set",1);
  /* check all path directories */
  pvec = str_split(path,':',&n);
  for (i=0,found=0; !found && i<n; i++) {
    p = str_file_path_new(pvec[i],prog);
    if (access(p,F_OK)==0 && access(p,X_OK)==0) {
      strcpy(prog,p);
      found = 1;
      }
    str_free(p);
    }
  vec_sfree(pvec,n);
  /* result check */
  if (!found)
    msg_error_f("cannot found executable file %s",1,prog);
  }

/* return string with path to the given executable program

   name - name of the program */
char* file_exec_path_find_new(char *name) {
  char *pv,**dr,*p,*path = NULL;
  unsigned i,n;
  short found;
  /* retrieve path variable */
  pv = getenv("PATH");
  if (!pv)
    msg_error("environment variable PATH is not set",1);
  /* check all path directories */
  dr = str_split(pv,':',&n);
  for (i=0,found=0; !found && i<n; i++) {
    p = str_file_path_new(dr[i],name);
    if (access(p,F_OK)==0 && access(p,X_OK)==0) {
      path = p;
      p = NULL;
      found = 1;
      }
    str_free(p);
    }
  vec_sfree(dr,n);
  /* result check */
  if (!found)
    msg_error_f("cannot found executable file %s",1,name);
  return(path);
  }

/* find path to given executable file

   name - name of the executable file */
void file_exec_path(char *name) {
  char *dir,*path;
  switch (file_exec_path_type(name)) {
    case PATH_ABSOLUTE:
      file_exec_path_check(name);
      break;
    case PATH_RELATIVE:
      dir = dir_current_new();
      path = str_file_path_new(dir,name);
      file_exec_path_check(path);
      strcpy(name,path);
      str_free(path);
      str_free(dir);
      break;
    case PATH_NAMEONLY:
      file_exec_path_find(name);
      break;
    }
  }

/* return string path to given executable file

   name - name of the executable file */
char* file_exec_path_new(char *name) {
  char *dir,*path = NULL;
  switch (file_exec_path_type(name)) {
    case PATH_ABSOLUTE:
      file_exec_path_check(name);
      path = str_copy_new(name);
      break;
    case PATH_RELATIVE:
      dir = dir_current_new();
      path = str_file_path_new(dir,name);
      file_exec_path_check(path);
      str_free(dir);
      break;
    case PATH_NAMEONLY:
      path = file_exec_path_find_new(name);
      break;
    }
  return(path);
  }

/* -------------------------------------------------------------------------- */
