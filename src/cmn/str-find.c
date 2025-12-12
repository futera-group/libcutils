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
#include <string.h>
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* find word in a file, stop at first occurence

   file - pointer to open file 
   word - the word which is searching for
   line - line which includes the word */
int str_ffind(FILE *file, char *word, char *line) {
  char *s;
  for (s=str_read_line_new(file); s;
       s=str_free(s),s=str_read_line_new(file)) {
    if (strstr(s,word)) {
      strcpy(line,s);
      str_free(s);
      return(1);
      }
    }
  return(0);
  }

/* find word in a file, stop at first occurence and return the line

   file - pointer to open file 
   word - the word which is searching for */
char* str_ffind_new(FILE *file, char *word) {
  char *s;
  for (s=str_read_line_new(file); s;
       s=str_free(s),s=str_read_line_new(file)) {
    if (strstr(s,word))
      return(s);
    }
  return(NULL);
  }

/* find word in a file, at the beggining of line, stop at first occurence

   file - pointer to open file 
   word - the word which is searching for
   line - line which includes the word */
int str_ffind_b(FILE *file, char *word, char *line) {
  char *s;
  for (s=str_read_line_new(file); s;
       s=str_free(s),s=str_read_line_new(file)) {
    if (str_sub_bfind(s,word)) {
      strcpy(line,s);
      str_free(s);
      return(1);
      }
    }
  return(0);
  }

/* find word in a file, at the beggining of line, return the line at first occ.

   file - pointer to open file 
   word - the word which is searching for */
char* str_ffind_b_new(FILE *file, char *word) {
  char *s;
  for (s=str_read_line_new(file); s;
       s=str_free(s),s=str_read_line_new(file)) {
    if (str_sub_bfind(s,word))
      return(s);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
