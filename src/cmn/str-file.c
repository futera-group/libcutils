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
#include <stdlib.h>
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* cut name of file from the path without extension
   
   path - path including file name (input)
   base - name of the file (output) */
int str_file_base(char *path, char *base) {
  int i,j,k;
  unsigned n;
  /* locate the base name */
  n = str_length(path);
  for (i=n-1; i>0; i--)
    if (path[i]=='/')
      break;
  for (j=n-1; j>0; j--)
    if (path[j]=='.')
      break;
  /* copy the base name */
  base[0] = '\0';
  if (j<=i)
    sprintf(base,"%s",path+i+1);
  else if (j>(i+1)) {
    if (i) {
      for (k=(i+1); k<j; k++)
        base[k-i-1] = path[k];
      base[j-i-1] = '\0';
      }
    else {
      for (k=0; k<j; k++)
        base[k] = path[k];
      base[j] = '\0';
      }
    }
  if (base[0]=='\0')
    return(0);
  return(1);
  }

/* return base nae of a file
   
   path - path including file name (input)
   name - name of the file (output) */
char* str_file_base_new(char *path) {
  int i,j,k;
  char *base = NULL;
  unsigned n;
  /* locate the base name */
  n = str_length(path);
  for (i=n-1; i>0; i--)
    if (path[i]=='/')
      break;
  for (j=n-1; j>0; j--)
    if (path[j]=='.')
      break;
  /* copy the base name */
  if (j<=i) {
    base = str_copy_new(path+i+1);
    }
  else if (j>(i+1)) {
    if (i) {
      base = str_new(j-i-1);
      for (k=(i+1); k<j; k++)
        base[k-i-1] = path[k];
      base[j-i-1] = '\0';
      }
    else {
      base = str_new(j);
      for (k=0; k<j; k++)
        base[k] = path[k];
      base[j] = '\0';
      }
    }
  return(base);
  }

/* -------------------------------------------------------------------------- */

/* cut name of file from the path
   
   path - path including file name (input)
   name - name of the file (output) */
int str_file_name(char *path, char *name) {
  int i;
  unsigned n;
  n = str_length(path);
  for (i=n-1; i>0; i--)
    if (path[i]=='/')
      break;
  if (i<0)
    name[0] = '\0';
  else {
    if (i==0) {
      if (path[0]=='/')
        sprintf(name,"%s",path+1);
      else
        sprintf(name,"%s",path);
      }
    else
      sprintf(name,"%s",path+i+1);
    }
  if (name[0]=='\0')
    return(0);
  return(1);
  }

/* cut name of file from the path and return it
   
   path - path including file name */
char* str_file_name_new(char *path) {
  short found;
  unsigned n;
  int i;
  char *file = NULL;
  if (path) {
    /* length of the string */
    n = str_length(path);
    for (i=n-1,found=0; i>0; i--)
      if (path[i]=='/') {
        found = 1;
        break;
        }
    /* copy directory name */
    if (!found && path[0]) 
      file = str_copy_new(path);
    else if (found)
      file = str_copy_new(path+i+1);
    }
  return(file);
  }

/* -------------------------------------------------------------------------- */

/* cut file extension from the path
   
   path - path including file name (input)
   ext  - file extension (output) */
int str_file_ext(char *path, char *ext) {
  int i,j;
  unsigned n;
  /* locate the base name */
  n = str_length(path);
  for (i=n-1; i>0; i--)
    if (path[i]=='/')
      break;
  for (j=n-1; j>0; j--)
    if (path[j]=='.')
      break;
  /* copy the extension */
  ext[0] = '\0';
  if (j<=i)
    return(0);
  sprintf(ext,"%s",path+j+1);
  return(1);
  }

/* return extension of a file
   
   path - path including file name */
char* str_file_ext_new(char *path) {
  int i,j;
  unsigned n;
  char *ext = NULL;
  /* locate the base name */
  n = str_length(path);
  for (i=n-1; i>0; i--)
    if (path[i]=='/')
      break;
  for (j=n-1; j>0; j--)
    if (path[j]=='.')
      break;
  /* copy the extension */
  if (j>i)
    ext = str_copy_new(path+j+1);
  return(ext);
  }

/* -------------------------------------------------------------------------- */

/* cut name of directory from the path without extension
   
   path - path including file name (input)
   name - name of the file (output) */
int str_file_dir(char *path, char *name) {
  int i,j;
  unsigned n;
  name[0]='\0';
  n=str_length(path);
  for (i=n-1; i>0; i--)
    if (path[i]=='/')
      break;
  for (j=0; j<i; j++)
    name[j]=path[j];
  name[j]='\0';
  if (name[0]=='\0')
    return(0);
  return(1);
  }

/* cut name of directory from the path without extension and return it
   
   path - path including file name */
char* str_file_dir_new(char *path) {
  int i,j;
  unsigned n;
  char *dir = NULL;
  if (path) {
    /* length of the string */
    n = str_length(path);
    for (i=n-1; i>0; i--)
      if (path[i]=='/')
        break;
    /* copy directory name */
    if (i>0) {
      dir = str_new(i);
      for (j=0; j<i; j++)
        dir[j] = path[j];
      dir[i]='\0';
      }
    }
  return(dir);
  }

/* -------------------------------------------------------------------------- */

/* create full path from directory and file name

   dir  - directory path
   name - name of the file
   path - storage for the path */
int str_file_path(char *dir, char *name, char *path) {
  if (dir && name)
    sprintf(path,"%s/%s",dir,name);
  else if (!dir && name)
    sprintf(path,"%s",name);
  else 
    return(0);
  return(1);
  }

/* return full path to a file

   dir  - directory path
   name - name of the file */
char *str_file_path_new(char *dir, char *name) {
  unsigned i,n1,n2;
  char *path = NULL;
  if (dir && name) {
    /* string lengths */
    n1 = str_length(dir);
    n2 = str_length(name);
    /* allocate memory */
    path = str_new(n1+n2+1);
    /* concatenated string */
    for (i=0; i<n1; i++)
        path[i] = dir[i];
      path[n1] = '/';
      for (i=0; i<n2; i++)
        path[i+n1+1] = name[i];
      path[n1+n2+1] = '\0';
    }
  else if (!dir && name)
    path = str_copy_new(name);
  return(path);
  }

/* -------------------------------------------------------------------------- */
