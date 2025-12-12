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

/* print text to the center of the line 
   
   str - string which should be printed */
void print_center(char *str) {
  unsigned int n,slen;
  slen = str_length(str);
  if (slen>=PR_LINE_WIDTH)
    printf("%s\n",str);
  else {
    n = (PR_LINE_WIDTH-slen)/2+slen;
    printf("%*s\n",n,str);
    }
  }

/* print text to the center of the line to file
   
   file - pointer to open file
   str  - string which should be printed */
void print_fcenter(FILE *file, char *str) {
  unsigned int n,slen;
  slen = str_length(str);
  if (slen>=PR_LINE_WIDTH)
    fprintf(file,"%s\n",str);
  else
    {
    n = (PR_LINE_WIDTH-slen)/2+slen;
    fprintf(file,"%*s\n",n,str);
    }
  }

/* print text to the center of the line and fill spaces with given char
   
   str - string which should be printed
   c   - the char for filling */
void print_ccenter(char *str, char c) {
  unsigned i,n,slen;
  slen = str_length(str);
  if (slen>=PR_LINE_WIDTH)
    printf("%s\n",str);
  else {
    n = (PR_LINE_WIDTH-slen)/2-1;
    for (i=0; i<n; i++)
      printf("%c",c);
    printf(" %s ",str);
    n = PR_LINE_WIDTH-n-slen-2;
    for (i=0; i<n; i++)
      printf("%c",c);
    printf("\n");
    }
  }

/* print text to the center of the line and fill spaces with given char
   
   file - pointer to open file
   str  - string which should be printed
   c    - the char for filling */
void print_fccenter(FILE *file, char *str, char c) {
  unsigned i,n,slen;
  slen = str_length(str);
  if (slen>=PR_LINE_WIDTH)
    fprintf(file,"%s\n",str);
  else {
    n = (PR_LINE_WIDTH-slen)/2-1;
    for (i=0; i<n; i++)
      fprintf(file,"%c",c);
    fprintf(file," %s ",str);
    n = PR_LINE_WIDTH-n-slen-2;
    for (i=0; i<n; i++)
      fprintf(file,"%c",c);
    fprintf(file,"\n");
    }
  }

/* -------------------------------------------------------------------------- */
