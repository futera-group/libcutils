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

#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* count number of specific characters in given string
 
   s - the string
   x - the character */
unsigned str_char_count(char *s, char x) {
  unsigned i,n;
  n = 0;
  for (i=0; i<str_length(s); i++)
    if (s[i]==x)
      n++;
  return(n);
  }

/* -------------------------------------------------------------------------- */

/* get the first blank character in the string
 
   s - pointer to the string
   n - id of the first character (output) */
char str_char_first_b(char *s, unsigned *n) {
  long int i = 0;
  char c='\0';
  if (n) 
    (*n) = 0;
  while (i<str_length(s)) {
    if (s[i]==' ' || s[i]=='\n' || s[i]=='\t' || s[i]=='\0')
      break;
    i++;
    }
  if (i<str_length(s)) {
    c = s[i];
    if (n) 
      (*n) = i;
    }
  return(c);
  }

/* Get the first non-blank character in the string 
 
   s - the input string
   n - ID of the first non-blank character (output,optional) */
char str_char_first_nb(char *s, unsigned *n) {
  long int i = 0;
  char c = '\0';
  if (n)
    (*n) = 0;
  while (i<str_length(s)) {
    if (s[i]!=' ' && s[i]!='\n' && s[i]!='\t' && s[i]!='\0')
      break;
    i++;
    }
  if (i<str_length(s)) {
    c = s[i];
    if (n)
      (*n) = i;
    }
  return(c);
  }

/* get the last non-blank character in the string
 
   s - pointer to the string
   n - id of the last character (output) */
char str_char_last_nb(char *s, unsigned *n) {
  long int i;
  char c='\0';
  if (n) 
    (*n) = 0;
  if ((i=str_length(s))) {
    while (i>=0) {
      if (s[i]!=' ' && s[i]!='\n' && s[i]!='\t' && s[i]!='\0')
        break;
      i--;
      }
    if (i>=0) {
      c = s[i];
      if (n) 
        (*n) = i;
      }
    }
  return(c);
  }

/* -------------------------------------------------------------------------- */
