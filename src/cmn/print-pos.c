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
#include "cmn/print.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* print text to given position in the line and fill spaces with given char
   
   str - string which should be printed
   c   - the char for filling */
void print_cpos(char *str, unsigned p, char c) {
  unsigned i,slen;
  slen = str_length(str);
  if (p<2)
    printf("%s",str);
  else {
    for (i=0; i<(p-1); i++)
      printf("%c",c);
    printf(" %s",str);
    }
  if ((slen+p+2)<=PR_LINE_WIDTH) {
    printf(" ");
    for (i=0; i<(PR_LINE_WIDTH-slen-p-1); i++)
      printf("%c",c);
    }
  printf("\n");
  }

/* -------------------------------------------------------------------------- */
