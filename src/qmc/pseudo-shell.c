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

#include <cmn/message.h>
#include <cmn/string.h>
#include "qmc/pseudo.h"

/* -------------------------------------------------------------------------- */

/* return ID of electronic shell
 
   s - name of the electronic shell */
short pseudo_shell_id(char *s) {
  str_trim(s);
  if (s[0]=='S')
    return(PSEUDO_L_S);
  if (s[0]=='P')
    return(PSEUDO_L_P);
  if (s[0]=='D')
    return(PSEUDO_L_D);
  if (s[0]=='F')
    return(PSEUDO_L_F);
  if (s[0]=='G')
    return(PSEUDO_L_G);
  else
    msg_error("invalid symbol of electronic shell",1);
  return(0);
  }

/* return name of electronic shell
 
   c - code of the electronic shell */
char *pseudo_shell_name(short c) {
  static char name[10];
  switch (c) {
    case PSEUDO_L_S: sprintf(name,"S"); break;
    case PSEUDO_L_P: sprintf(name,"P"); break;
    case PSEUDO_L_D: sprintf(name,"D"); break;
    case PSEUDO_L_F: sprintf(name,"F"); break;
    case PSEUDO_L_G: sprintf(name,"G"); break;
    default: name[0]='\0'; break;
    }
  return(name);
  }

/* -------------------------------------------------------------------------- */
