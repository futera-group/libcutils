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
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* convert string to lowcase chars, the string is changed

   str - the string which will be changed */
char* str_lowcase(char *str) {
  unsigned i;
  for (i=0; i<str_length(str); i++)
    if (str[i]>='A' && str[i]<='Z')
      str[i] = str[i]-'A'+'a';
  return(str);
  }

/* convert string to lowcase chars, the string is not changed

   s_in  - input string
   s_out - output string */
void str_lowcase_copy(char *s_in, char *s_out) {
  unsigned i;
  for (i=0; i<str_length(s_in); i++)
    if (s_in[i]>='A' && s_in[i]<='Z')
      s_out[i] = s_in[i]-'A'+'a';
    else
      s_out[i] = s_in[i];
  s_out[i] = '\0';
  }

/* create new lowercase version of input string

   s - the input string */
char* str_lowcase_new(char *s) {
  char *str = NULL;
  str = str_copy_new(s);
  str_lowcase(str);
  return(str);
  }

/* -------------------------------------------------------------------------- */

/* convert string to lowcase chars instead of the first which is capital

   str  - input string */
char* str_Lowcase(char *str) {
  unsigned i,slen;
  slen = str_length(str);
  if (slen==0)
    return(str);
  if (str[0]>='a' && str[0]<='z')
    str[0] = str[0]-'a'+'A';
  for (i=1; i<slen; i++)
    if (str[i]>='A' && str[i]<='Z')
      str[i] = str[i]-'A'+'a';
  return(str);
  }

/* convert string to lowcase chars instead of the first which is capital

   s_in  - input string
   s_out - output string */
void str_Lowcase_copy(char *s_in, char *s_out) {
  unsigned i,slen;
  slen=str_length(s_in);
  if (slen==0) return;
  if (s_in[0]>='a' && s_in[0]<='z')
    s_out[0] = s_in[0]-'a'+'A';
  else
    s_out[0] = s_in[0];
  for (i=1; i<slen; i++)
    if (s_in[i]>='A' && s_in[i]<='Z')
      s_out[i] = s_in[i]-'A'+'a';
    else
      s_out[i] = s_in[i];
  s_out[i] = '\0';
  }

/* create new lowercase version of input string with first capital

   s - the input string */
char* str_Lowcase_new(char *s) {
  char *str = NULL;
  str = str_copy_new(s);
  str_Lowcase(str);
  return(str);
  }

/* -------------------------------------------------------------------------- */

/* convert string to upcase chars, the string is changed

   str - the string which will be changed */
char* str_upcase(char *str) {
  unsigned i;
  for (i=0; i<str_length(str); i++)
    if (str[i]>='a' && str[i]<='z')
      str[i] = str[i]-'a'+'A';
  return(str);
  }


/* convert string to upcase chars, the string is not changed

   s_in  - input string
   s_out - output string */
void str_upcase_copy(char *s_in, char *s_out) {
  unsigned i;
  for (i=0; i<str_length(s_in); i++)
    if (s_in[i]>='a' && s_in[i]<='z')
      s_out[i] = s_in[i]-'a'+'A';
    else 
      s_out[i] = s_in[i];
  s_out[i] = '\0';
  }

/* create new upper-case version of input string

   s - the input string */
char* str_upcase_new(char *s) {
  char *str = NULL;
  str = str_copy_new(s);
  str_upcase(str);
  return(str);
  }

/* -------------------------------------------------------------------------- */
