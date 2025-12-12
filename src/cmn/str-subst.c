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

/* substitute all specified characters in the string
 
   s  - the string
   c1 - the character that is to be substituted
   c2 - the new character */
void str_subst(char *s, char c1, char c2) {
  unsigned i;
  for (i=0; i<str_length(s); i++)
    if (s[i]==c1)
      s[i] = c2;
  }

/* substitute first specified characters in the string
 
   s  - the string
   c1 - the character that is to be substituted
   c2 - the new character */
void str_subst_one(char *s, char c1, char c2) {
  unsigned i;
  for (i=0; i<str_length(s); i++)
    if (s[i]==c1) {
      s[i] = c2;
      break;
    }
  }

/* substitute first specified characters in the string by u-int
 
   s - the string
   c - the character that is to be substituted
   n - the number value */
void str_subst_one_ui(char **s, char c, unsigned n) {
  char *s1,*s2;
  unsigned i,j,ns,nlen;
  if (s && *s && (*s)[0]) {
    /* initialization */
    ns = str_length(*s);
    s1 = str_new(ns);
    s2 = str_new(ns);
    /* first part of the string */
    for (i=0; i<ns; i++) {
      s1[i] = (*s)[i];
      if ((*s)[i] == c)
        break;
      }
    s1[i] = '\0';
    /* second part of the string */
    j = 0;
    while (++i < ns)
      s2[j++] = (*s)[i];
    s2[j] = '\0';
    /* memory reallocation */
    nlen = str_unum_len(n);
    if (nlen > 1)
      *s = (char*)realloc(*s,ns+nlen);
    /* substitution */
    sprintf(*s,"%s%d%s",s1,n,s2);
    (*s)[ns+nlen-1] = '\0';
    /* clean memory */
    str_free(s1);
    str_free(s2);
    }
  }

/* -------------------------------------------------------------------------- */
