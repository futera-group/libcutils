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

/* compare two strings and return 0 or 1 if they are different or not

  s1,s2 - the strings */
int str_compare(char *s1, char *s2) {
  int i,l1,l2;
  l1 = str_length(s1);
  while (l1>0 && (s1[l1-1]==' ' || s1[l1-1]=='\t' || s1[l1-1]=='\n'))
    l1--;
  l2 = str_length(s2);
  while (l2>0 && (s2[l2-1]==' ' || s2[l2-1]=='\t' || s2[l2-1]=='\n'))
    l2--;
  if (l1==0 && l2==0)
    return(1);
  if (l1!=l2)
    return(0);
  i = 0;
  while (i<l1) {
    if (s1[i]!=s2[i])
      return(0);
    i++;
    }
  return(1);
  }

/* case-insensitive comparison of two strings 

  s1,s2 - the strings */
int str_compare_nc(char *s1, char *s2) {
  int i,l1,l2;
  char c1,c2;
  l1 = str_length(s1);
  while (l1>0 && (s1[l1-1]==' ' || s1[l1-1]=='\t' || s1[l1-1]=='\n'))
    l1--;
  l2 = str_length(s2);
  while (l2>0 && (s2[l2-1]==' ' || s2[l2-1]=='\t' || s2[l2-1]=='\n'))
    l2--;
  if (l1==0 && l2==0)
    return(1);
  if (l1!=l2)
    return(0);
  i = 0;
  while (i<l1) {
    c1 = (s1[i]>='A' && s1[i]<='Z' ? s1[i]-'A'+'a' : s1[i]);
    c2 = (s2[i]>='A' && s2[i]<='Z' ? s2[i]-'A'+'a' : s2[i]);
    if (c1!=c2)
      return(0);
    i++;
    }
  return(1);
  }

/* -------------------------------------------------------------------------- */
