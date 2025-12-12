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

/* print text and int number with fix line width 

   str  - string which should be printed
   data - number which should be printed
   num  - total length of printed text */
void print_iwfix(char *str, int data, unsigned num) {
  unsigned int slen,len;
  slen = str_length(str);
  len = (num-slen > 0 ? num-slen : 0);
  printf("%s%*d\n",str,len,data);
  }

/* print text and int number with fix line width to file

   file - pointer to open file
   str  - string which should be printed
   data - number which should be printed
   num  - total length of printed text */
void print_fiwfix(FILE *file, char *str, int data, unsigned num) {
  unsigned int slen,len;
  slen = str_length(str);
  len = (num-slen > 0 ? num-slen : 0);
  fprintf(file,"%s%*d\n",str,len,data);
  }

/* -------------------------------------------------------------------------- */

/* print text and unsigned int number with fix line width 

   str  - string which should be printed
   data - number which should be printed
   num  - total length of printed text */
void print_uwfix(char *str, unsigned data, unsigned num) {
  unsigned int slen,len;
  slen = str_length(str);
  len = (num-slen > 0 ? num-slen : 0);
  printf("%s%*u\n",str,len,data);
  }

/* print text and unsigned int number with fix line width to file

   file - pointer to open file
   str  - string which should be printed
   data - number which should be printed
   num  - total length of printed text */
void print_fuwfix(FILE *file, char *str, unsigned data, unsigned num) {
  unsigned int slen,len;
  slen = str_length(str);
  len = (num-slen > 0 ? num-slen : 0);
  fprintf(file,"%s%*u\n",str,len,data);
  }

/* -------------------------------------------------------------------------- */

/* print text and unsigned long int number with fix line width 

   str  - string which should be printed
   data - number which should be printed
   num  - total length of printed text */
void print_luwfix(char *str, unsigned long data, unsigned num) {
  unsigned int slen,len;
  slen = str_length(str);
  len = (num-slen > 0 ? num-slen : 0);
  printf("%s%*lu\n",str,len,data);
  }

/* print text and unsigned long int number with fix line width to file

   file - pointer to open file
   str  - string which should be printed
   data - number which should be printed
   num  - total length of printed text */
void print_fluwfix(FILE *file, char *str, unsigned long data, unsigned num) {
  unsigned int slen,len;
  slen = str_length(str);
  len = (num-slen > 0 ? num-slen : 0);
  fprintf(file,"%s%*lu\n",str,len,data);
  }

/* -------------------------------------------------------------------------- */

/* print text and unsigned int in octal format with fix line width 

   str  - string which should be printed
   data - number which should be printed
   num  - total length of printed text */
void print_lowfix(char *str, unsigned long data, unsigned num) {
  unsigned int slen,len;
  slen = str_length(str);
  len = (num-slen > 0 ? num-slen : 0);
  printf("%s%*lo\n",str,len,data);
  }

/* print text and unsigned int in octal format with fix line width to file

   file - pointer to open file
   str  - string which should be printed
   data - number which should be printed
   num  - total length of printed text */
void print_flowfix(FILE *file, char *str, unsigned long data, unsigned num) {
  unsigned int slen,len;
  slen = str_length(str);
  len = (num-slen > 0 ? num-slen : 0);
  fprintf(file,"%s%*lo\n",str,len,data);
  }

/* -------------------------------------------------------------------------- */

/* print text and real number with fix line width

   str  - string which should be printed
   frmt - format for real number printing
   data - the real number 
   num  - total length of printed text */
void print_fwfix(char *str, char *frmt, double data, unsigned num) {
  char snum[PR_LINE_WIDTH];
  unsigned slen,len;
  sprintf(snum,frmt,data);
  slen = str_length(str);
  len = (num-slen > 0 ? num-slen : 0);
  printf("%s%*s\n",str,len,snum);
  }

/* print text and real number with fix line width to file

   file - pointer to open file
   str  - string which should be printed
   frmt - format for real number printing
   data - the real number 
   num  - total length of printed text */
void print_ffwfix(FILE *file, char *str, char *frmt, double data,
  unsigned num) {
  char snum[PR_LINE_WIDTH];
  unsigned slen,len;
  sprintf(snum,frmt,data);
  slen = str_length(snum);
  len = (num-slen > 0 ? num-slen : 0);
  printf("%s%*s\n",str,len,snum);
  }

/* -------------------------------------------------------------------------- */

/* print text and string with fix line width 

   str  - string which should be printed
   data - text which should be printed after str
   num  - total length of printed text */
void print_swfix(char *str, char *data, unsigned num) {
  unsigned int slen,len;
  slen = str_length(str);
  len = (num-slen > 0 ? num-slen : 0);
  printf("%s%*s\n",str,len,data);
  }

/* print text and string with fix line width to file

   file - pointer to open file
   str  - string which should be printed
   data - text which should be printed after str
   num  - total length of printed text */
void print_fswfix(FILE *file, char *str, char *data, unsigned num) {
  unsigned int slen,len;
  slen = str_length(str);
  len = (num-slen > 0 ? num-slen : 0);
  fprintf(file,"%s%*s\n",str,len,data);
  }

/* -------------------------------------------------------------------------- */
