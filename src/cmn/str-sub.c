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

/* copy substring from specific interval

   source - string from which is copied
   in     - string to which is copied
   n      - first character of the substring
   m      - last character of the substring */
void str_sub_copy(char *source, char *in, unsigned n, unsigned m) {
  unsigned i,j=0;
  if (n<m) {
    for (i=n; i<m; i++) {
      in[j] = source[i];
      j++;
      }
    in[j] = '\0';
    }
  }

/* copy substring after specified deliminer
 
   s - string from which is copied
   d - string to which is copie
   c - deliminer character */
void str_sub_copy_adel(char *s, char *d, char c) {
  unsigned i=0,j=0,n;
  n = str_length(s);
  while (i<n && s[i]!=c)
    i++;
  i++;
  while (i<n) {
    d[j] = s[i];
    i++;
    j++;
    }
  d[j] = '\0';
  }

/* copy substring before specified deliminer
 
   s - string from which is copied
   d - string to which is copie
   c - deliminer character */
void str_sub_copy_bdel(char *s, char *d, char c) {
  unsigned i=0,n;
  n = str_length(s);
  while (i<n && s[i]!=c) {
    d[i] = s[i];
    i++;
    }
  d[i] = '\0';
  }

/* -------------------------------------------------------------------------- */

/* find substring at the begin of string (return 1 if found)

   where - string where to search
   what  - substring which is searched for */
int str_sub_bfind(char *where ,char *what) {
  unsigned i=0,n,m;
  n = str_length(where);
  m = str_length(what);
  n = (n<m ? n : m);
  while (i<n) {
    if (where[i]!=what[i])
      return(0);
    i++;
    }
  return(1);
  }

/* -------------------------------------------------------------------------- */
