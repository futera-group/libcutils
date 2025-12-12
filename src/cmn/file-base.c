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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "cmn/message.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* open file and return its pointer

   file_name - name of the file which should be opened
   mode      - opening mode (e.i. "r","w","rw",...) */
FILE *file_open(char *file_name, char *mode) {
  FILE *file = NULL;
  file = fopen(file_name,mode);
  if (!file)
    msg_error_f("cannot open file \"%s\"",1,file_name);
  return(file);
  }

/* open temporary file */
FILE *file_open_tmp(void) {
  FILE *file = NULL;
  file = tmpfile();
  if (!file)
    msg_error("cannot open temporary file",1);
  return(file);
  }

/* associate open file stream with different file
 
   name - name of the new file
   mode - opening mode ("r","w","rw",...)
   fs   - the file stream */
FILE *file_reopen(char *name, char *mode, FILE *fs) {
  FILE *f = NULL;
  f = freopen(name,mode,fs);
  if (!f)
    msg_error_f("cannot reopen stream with file \"%s\"",1,name);
  return(f);
  }

/* close opened file

   file - the file which should be closed */
void file_close(FILE *file) {
  if (fclose(file)!=0)
    msg_warn("file closing failed");
  }

/* -------------------------------------------------------------------------- */

/* get information about file

   file  - string file name
   fstat - pointer to stat structure 
   code  - exit if not 0 when an error occurs */
int file_stat(char *file, struct stat *fstat, short code) {
  char error[80];
  if (lstat(file,fstat)<0) {
    switch (errno) {
      case EBADF:        sprintf(error,"bad file descriptor");        break;
      case ENOENT:       sprintf(error,"no such file or directory");  break;
      case ENOTDIR:      sprintf(error,"not a directory");            break;
      case ELOOP:        sprintf(error,"too many symbolic links");    break;
      case EFAULT:       sprintf(error,"bad address");                break;
      case EACCES:       sprintf(error,"permission denied");          break;
      case ENOMEM:       sprintf(error,"insufficient kernel memory"); break;
      case ENAMETOOLONG: sprintf(error,"too long name");              break;
      default:           sprintf(error,"unknown error");              break;
      }
    msg_error_f("cannot get info about %s (%s)",1,file,error);
    return(0);
    }
  return(1);
  }

/* detect if given file exists

   path - path to the file */
int file_exist(char *path) {
  struct stat st;
  if (lstat(path,&st)<0 && errno==ENOENT)
    return(0);
  return(1); 
  }

/* -------------------------------------------------------------------------- */

/* create symbolic link to file saved in the filesystem
 
   f1 - the original file name
   f2 - the symbolic link name */
void file_link(char *f1, char *f2) {
  char msg[80];
  if (symlink(f1,f2)==-1) {
    switch (errno) {
      case EPERM:
      case EACCES:       sprintf(msg,"permission denied");               break;
      case EDQUOT:       sprintf(msg,"disk quota exceeded");             break;
      case EEXIST:       sprintf(msg,"file already exists");             break;
      case EFAULT:       sprintf(msg,"out of address space");            break;
      case EIO:          sprintf(msg,"I/O error");                       break;
      case ELOOP:        sprintf(msg,"too many symlinks in path");       break;
      case ENAMETOOLONG: sprintf(msg,"too long name");                   break;
      case ENOENT:       sprintf(msg,"no such file");                    break;
      case ENOSPC:       sprintf(msg,"no space");                        break;
      case ENOTDIR:      sprintf(msg,"part of path is not a directory"); break;
      case EROFS:        sprintf(msg,"read-only filesystem");            break;
      }
    msg_error_f("cannot create link to file %s (%s)",1,f1,msg);
    }
  }

/* move or rename file saved in the specified path
 
   f1 - the original file name
   f2 - the new file name */
void file_move(char *f1, char *f2) {
  char msg[80];
  if (rename(f1,f2)==-1) {
    switch (errno) {
      case EPERM:
      case EACCES:       sprintf(msg,"permission denied");               break;
      case EDQUOT:       sprintf(msg,"disk quota exceeded");             break;
      case EFAULT:       sprintf(msg,"out of address space");            break;
      case EINVAL:       sprintf(msg,"invalid path");                    break;
      case EISDIR:       sprintf(msg,"is a directory");                  break;
      case EIO:          sprintf(msg,"I/O error");                       break;
      case ELOOP:        sprintf(msg,"too many symlinks in path");       break;
      case ENAMETOOLONG: sprintf(msg,"too long name");                   break;
      case ENOENT:       sprintf(msg,"no such file");                    break;
      case ENOSPC:       sprintf(msg,"no space");                        break;
      case ENOTDIR:      sprintf(msg,"part of path is not a directory"); break;
      case ENOTEMPTY:    sprintf(msg,"directory is not empty");          break;
      case EROFS:        sprintf(msg,"read-only filesystem");            break;
      case EXDEV:        sprintf(msg,"different logical devices");       break;
      }
    msg_error_f("cannot move/rename file %s (%s)",1,f1,msg);
    }
  }

/* copy one file to another
  
   p1 - path to the source file
   p2 - path to the destination file */
void file_copy(char *p1, char *p2) {
  char *line;
  FILE *f1,*f2; 
  /* open the files */
  f1 = file_open(p1,"r");
  f2 = file_open(p2,"w");
  /* copy contents of the file */
  for (line=str_read_line_new(f1); line; 
       line=str_free(line), line=str_read_line_new(f1)) 
    fprintf(f2,"%s",line);
  str_free(line);
  /* close the files */
  file_close(f1);
  file_close(f2);
  }

/* merge two files to new one

   p1 - path to the first source file
   p2 - path to the second source file
   p3 - path to the desctination file */
void file_merge(char *p1, char *p2, char *p3) {
  char *line;
  FILE *f1,*f2;
  /* destination file */
  f2 = file_open(p3,"w");
  /* copy contents of the first file */
  f1 = file_open(p1,"r");
  for (line=str_read_line_new(f1); line; 
       line=str_free(line), line=str_read_line_new(f1)) 
    fprintf(f2,"%s",line);
  str_free(line);
  file_close(f1);
  /* add contents of the second file */
  f1 = file_open(p2,"r");
  for (line=str_read_line_new(f1); line; 
       line=str_free(line), line=str_read_line_new(f1)) 
    fprintf(f2,"%s",line);
  str_free(line);
  file_close(f1);
  /* close the file */
  file_close(f2);
  }

/* remove a file

   path - path to the file */
void file_remove(char *path) {
  char msg[80];
  if (unlink(path)==-1) {
    switch (errno) {
      case EPERM:
      case EACCES:       sprintf(msg,"permission denied");               break;
      case EISDIR:       sprintf(msg,"is a directory");                  break;
      case EBUSY:        sprintf(msg,"file is busy");                    break;
      case EFAULT:       sprintf(msg,"out of address space");            break;
      case ENAMETOOLONG: sprintf(msg,"too long name");                   break;
      case ENOENT:       sprintf(msg,"no such file");                    break;
      case ENOTDIR:      sprintf(msg,"part of path is not a directory"); break;
      case ENOMEM:       sprintf(msg,"insufficient kernel memory");      break;
      case EROFS:        sprintf(msg,"read-only filesystem");            break;
      case ELOOP:        sprintf(msg,"too many symlinks in path");       break;
      case EIO:          sprintf(msg,"I/O error");                       break;
      }
    msg_error_f("cannot remove file %s (%s)",1,path,msg);
    }
  }

/* -------------------------------------------------------------------------- */

/* return size of specified file 
 
   name - name of the file */
long unsigned file_size(char *name) {
  struct stat fs;
  if (!file_stat(name,&fs,0))
    msg_error_f("cannot obtain size of \"%s\" file",1,name);
  return(fs.st_size);
  }

/* -------------------------------------------------------------------------- */
