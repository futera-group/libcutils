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
#include <cmn/message.h>
#include <cmn/string.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* identify quantum numbers of from the orbital type keyword
 
   key - orbital type keyword
   n   - principal quantum number (output)
   l   - angular quantum number (output) 
   m   - magnetic quantum number (output) */
void cp2k_orbital_type_id(char *key, int *n, int *l, int *m) {
  char s[80] = "\0";
  /* initialization */
  (*n) = (*l) = (*m) = 0;
  /* principal quantum number */
  if (sscanf(key,"%d%s",n,s)!=2)
    msg_error_f("invalid orbital type specification \"%s\"",1,key);
  if ((*n)<1)
    msg_error_f("invalid principal quantum number (%d)"
      " of \"%s\" orbital",1,*n,key);
  /* angular and magnetic quantum numbers */
  str_lowcase(s);
  switch (s[0]) {
    case 's': (*l) = 0; /* S fce */
      break;
    case 'p': (*l) = 1; /* P fce */
      if      (str_compare(s+1,"x")) (*m) = 0;
      else if (str_compare(s+1,"y")) (*m) = 1;
      else if (str_compare(s+1,"z")) (*m) = 2;
      else 
        msg_error_f("unknown P-type \"%s\" of orbital \"%s\"",1,key,s+1);
      break;
    case 'd': (*l) = 2; /* D fce */
      if      (str_compare(s+1,"x2")) (*m) = 0;
      else if (str_compare(s+1,"xy")) (*m) = 1;
      else if (str_compare(s+1,"xz")) (*m) = 2;
      else if (str_compare(s+1,"y2")) (*m) = 3;
      else if (str_compare(s+1,"yz")) (*m) = 4;
      else if (str_compare(s+1,"z2")) (*m) = 5;
      else 
        msg_error_f("unknown D-type \"%s\" of orbital \"%s\"",1,key,s+1);
      break;
    case 'f':  (*l) = 3; /* F fce */
      if      (str_compare(s+1,"x3"))  (*m) = 0;
      else if (str_compare(s+1,"x2y")) (*m) = 1;
      else if (str_compare(s+1,"x2z")) (*m) = 2;
      else if (str_compare(s+1,"xy2")) (*m) = 3;
      else if (str_compare(s+1,"xyz")) (*m) = 4;
      else if (str_compare(s+1,"xz2")) (*m) = 5;
      else if (str_compare(s+1,"y3"))  (*m) = 6;
      else if (str_compare(s+1,"y2z")) (*m) = 7;
      else if (str_compare(s+1,"yz2")) (*m) = 8;
      else if (str_compare(s+1,"z3"))  (*m) = 9;
      else 
        msg_error_f("unknown F-type \"%s\" of orbital \"%s\"",1,key,s+1);
      break;
    default:
      msg_error_f("unknown angular type \"%c\" of orbital \"%s\"",1,key,s[0]);
    }
  }

/* -------------------------------------------------------------------------- */
