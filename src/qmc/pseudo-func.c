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
#include "qmc/pseudo.h"

/* -------------------------------------------------------------------------- */

/* return ID of exchange-correlation functional
 
   s - name of the exchange-correlation functional */
short pseudo_func_id(char *s) {
  if (strstr(s,"Lee-Yang-Parr"))
    return(PSEUDO_XC_LYP);
  if (strstr(s,"Becke (1988)"))
    return(PSEUDO_XC_B88);
  return(0);
  }

/* return name of exchange-correlation functional
 
   c - code of the exchange-correlation functional */
char *pseudo_func_name(short c) {
  static char name[80]="\0";
  switch (c) {
    case PSEUDO_XC_LYP: sprintf(name,"Lee-Yang-Parr"); break;
    case PSEUDO_XC_B88: sprintf(name,"Becke (1988)"); break;
    default: name[0]='\0'; break;
    }
  return(name);
  }

/* -------------------------------------------------------------------------- */
