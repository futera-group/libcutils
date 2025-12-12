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

/* trim spaces from both sides of the string

   str - input string which will be changed */
void str_trim(char *str) {
  int i,first,last;
  unsigned slen;
  if (!str)
    return;
  slen = str_length(str);
  /* first non-blank character */
  first = 0;
  while (first<slen &&
        (str[first]==' ' || str[first]=='\t' || str[first]=='\n'))
    first++;
  /* last non-blank character */
  last = slen-1;
  while (last-1>=0 &&
        (str[last]==' ' || str[last]=='\t' || str[last]=='\n'))
    last--;
  /* trimming */
  if (first<=last) {
    for (i=first; i<=last; i++)
      str[i-first] = str[i];
    str[i-first] = '\0';
    }
  else
    str[0]='\0';
  }

/* trim spaces from both sides of the string, the string is not changed

   s_in - input string
   s_out - output string  */
void str_trim_copy(char *s_in, char *s_out) {
  int i,first,last;
  unsigned slen;
  slen = str_length(s_in);
  /* first non-blank character */
  first = 0;
  while (first<slen &&
        (s_in[first]==' ' || s_in[first]=='\t' || s_in[first]=='\n'))
    first++;
  /* last non-blank character */
  last = slen-1;
  while (last-1>=0 &&
        (s_in[last]==' ' || s_in[last]=='\t' || s_in[last]=='\n'))
    last--;
  /* trimming */
  if (first<=last) {
    for (i=first; i<=last; i++)
      s_out[i-first] = s_in[i];
    s_out[i-first] = '\0';
    }
  else
    s_out[0]='\0';
  }

/* trim string from both sides up to specified character

   str - input string which will be changed 
   mark - character used as a mark for trimming */
void str_trim_mark(char *str, char mark) {
  str_ltrim_mark(str,mark);
  str_rtrim_mark(str,mark);
  }

/* delete specific characters from the string

   si - string for processing (not modified)
   sc - string of characters to delete
   st - storage for the ouput modified string */
void str_trim_char(char *si, char *sc, char *st) {
  unsigned i,j,id,n;
  short found;
  id = 0;
  n = str_length(sc);
  for (i=0; i<str_length(si); i++) {
    found = 0;
    for (j=0; j<n; j++) {
      if (sc[j]==si[i]) {
        found = 1;
        break;
        }
      if (sc[j]=='\0')
        break;
      }
    if (!found)
      st[id++] = si[i];
    if (si[i]=='\0')
      break;
    }
  st[id] = '\0';
  }

/* -------------------------------------------------------------------------- */

/* trim spaces from right side of the string

   str - input string which will be changed */
void str_rtrim(char *str) {
  int last;
  unsigned slen;
  slen = str_length(str);
  /* last non-blank character */
  last = slen-1;
  while (last>0 && (str[last]==' ' || str[last]=='\t' || str[last]=='\n'))
    last--;
  /* trimming */
  if (last==0)
    str[0] = '\0';
  else if (last==(slen-1))
    str[last] = '\0';
  else
    str[last+1] = '\0';
  }

/* trim string from right side up to specified character

   str - input string which will be changed 
   mark - character used as a mark for trimming */
void str_rtrim_mark(char *str, char mark) {
  unsigned slen;
  int last;
  slen = str_length(str);
  /* find the marker */
  last = slen-1;
  while (last>0 && str[last]!=mark)
    last--;
  /* trimming */
  if (last>0) {
    if (last==(slen-1))
      str[last] = '\0';
    else
      str[last] = '\0';
    }
  else if (str[0]==mark)
    str[0] = '\0';
  }

/* -------------------------------------------------------------------------- */

/* trim spaces from left side of the string

   str - input string which will be changed */
void str_ltrim(char *str) {
  int i,first;
  unsigned slen;
  slen = str_length(str);
  /* first non-blank character */
  first = 0;
  while (first<slen &&
        (str[first]==' ' || str[first]=='\t' || str[first]=='\n'))
    first++;
  /* trimming */
  if (first<slen) {
    for (i=first; i<slen; i++)
      str[i-first] = str[i];
    str[i-first] = '\0';
    }
  else
    str[0]='\0';
  }

/* trim string from left side up to specified character

   str - input string which will be changed 
   mark - character used as a mark for trimming */
void str_ltrim_mark(char *str, char mark) {
  unsigned i,slen;
  int first=0;
  slen = str_length(str);
  while (first<slen && str[first++]!=mark);
  for (i=0; i<slen && first<slen; i++)
    str[i] = str[first++];
  str[i] = '\0';
  }

/* -------------------------------------------------------------------------- */
