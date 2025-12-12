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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cmn/message.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for a string

   n - length of the string */
char* str_new(unsigned n) {
  char *s = NULL;
  /* memory allocation */
  s = (char*)malloc((n+1)*sizeof(char));
  if (!s)
    msg_error("cannot allocate memory for a string",1);
  /* initialization */
  s[0] = '\0';
  return(s);
  }

/* free memory allocated for string
 
   s - the allocated string */
char* str_free(char *s) {
  if (s) {
    free(s);
    s = NULL;
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */

/* return lenght of a string 

   s - the string */
unsigned str_length(char *s) {
  unsigned n = 0;
  if (s)
    n = strlen(s);
  return(n);
  }

/* -------------------------------------------------------------------------- */

/* merge two strings

   s1,s2 - the strings to be merged
   t     - storage for the merged string */
void str_merge(char *s1, char *s2, char *t) {
  unsigned i,n1,n2;
  n1 = str_length(s1);
  n2 = str_length(s2);
  for (i=0; i<n1; i++)
    t[i] = s1[i];
  for (i=0; i<n2; i++)
    t[i+n1] = s2[i];
  t[n1+n2] = '\0';
  }

/* merge two strings to a new one 

   s1,s2 - the strings to be merged */
char* str_merge_new(char *s1, char *s2) {
  unsigned i,n1,n2;
  char *t = NULL;
  if (s1 && s2) {
    n1 = strlen(s1);
    n2 = strlen(s2);
    t = str_new(n1+n2);
    for (i=0; i<n1; i++)
      t[i] = s1[i];
    for (i=0; i<n2; i++)
      t[i+n1] = s2[i];
    t[n1+n2] = '\0';
    }
  else if (s1) 
    t = str_copy_new(s1);
  else if (s2)
    t = str_copy_new(s2);
  return(t);
  }

/* -------------------------------------------------------------------------- */

/* append string to dynamically allocated one

   s1 - the dynamically allocated string
   s2 - the appendix */
void str_append(char **s1, char *s2) {
  char *t = NULL;
  if (s1) {
    if ((*s1) && s2) {
      t = str_merge_new(*s1,s2);
      str_free(*s1);
      (*s1) = t;
      }
    else if (s2)
      (*s1) = str_copy_new(s2);
    }
  }

/* -------------------------------------------------------------------------- */

/* formatted print to new dynamically-allocated string
 
   s - string template with formatting flags */
char* str_sprintf(char *s, ...) {
  char t[2],*x;
  int n;
  va_list ap;
  va_start(ap,s);
  n = vsnprintf(t,1,s,ap);
  va_end(ap);
  x = str_new(n+1);
  va_start(ap,s);
  n = vsnprintf(x,n+1,s,ap);
  va_end(ap);
  return(x);
  }

/* -------------------------------------------------------------------------- */

/* copy a string
 
   s1 - storage for a copy 
   s2 - the original string */
void str_copy(char *s1, char *s2) {
  strcpy(s1,s2);
  }

/* creates new copy of the character string
 
   s - the string */
char* str_copy_new(char *s) {
  unsigned n;
  char *t = NULL;
  if (s) {
    n = strlen(s);
    t = str_new(n);
    strncpy(t,s,n);
    t[n] = '\0';
    }
  return(t);
  }

/* -------------------------------------------------------------------------- */
