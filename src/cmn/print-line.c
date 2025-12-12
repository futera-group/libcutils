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
#include "cmn/print.h"

/* -------------------------------------------------------------------------- */

/* print horizontal line */
void print_hline(void) {
  unsigned int i;
  for (i=0; i<PR_LINE_WIDTH; i++)
    printf("-");
  printf("\n");
  }

/* print horizontal line to file

   file - pointer to open file */
void print_fhline(FILE *file) {
  unsigned int i;
  for (i=0; i<PR_LINE_WIDTH; i++)
    fprintf(file,"-");
  fprintf(file,"\n");
  }

/* print horizontal line with specified chars and length 

   c - char for drawing line
   n - number of chars (lenght of the line) */
void print_Hline(char c, unsigned n) {
  unsigned int i;
  for (i=0; i<n; i++)
    printf("%c",c);
  printf("\n");
  }

/* print horizontal line with specified chars and length to file

   file - pointer to open file
   c    - char for drawing line
   n    - number of chars (lenght of the line) */
void print_fHline(FILE *file, char c, unsigned n) {
  unsigned int i;
  for (i=0; i<n; i++)
    fprintf(file,"%c",c);
  fprintf(file,"\n");
  }

/* -------------------------------------------------------------------------- */
