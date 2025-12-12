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

#ifndef ZF_LIB_CMN_PRINT_H
#define ZF_LIB_CMN_PRINT_H

#include <stdio.h>

/* -------------------------------------------------------------------------- */

/* symbolic constant */
#define PR_LINE_WIDTH 80

/* -------------------------------------------------------------------------- */

/* print horizontal line */
void print_hline(void);
/* print horizontal line to file */
void print_fhline(FILE*);
/* print horizontal line with specified chars and length */
void print_Hline(char, unsigned);
/* print horizontal line with specified chars and length to file*/
void print_fHline(FILE*, char, unsigned);

/* print text to given position in the line and fill spaces with given char */
void print_cpos(char*, unsigned, char);

/* print text to the center of the line */
void print_center(char*);
/* print text to the center of the line to file */
void print_fcenter(FILE*, char*);
/* print text to the center of the line and fill spaces with given char */
void print_ccenter(char*,char);
/* print text to the center of the line and fill spaces with given char */
void print_fccenter(FILE*, char*,char);

/* print text and int number to fix width */
void print_iwfix(char*, int, unsigned);
/* print text and int number to fix width to file */
void print_fiwfix(FILE*, char*, int, unsigned);

/* print text and unsigned int number to fix width */
void print_uwfix(char*, unsigned, unsigned);
/* print text and unsigned int number to fix width to file*/
void print_fuwfix(FILE*, char*, unsigned, unsigned);

/* print text and unsigned long int number to fix width */
void print_luwfix(char*, unsigned long, unsigned);
/* print text and unsigned long int number to fix width to file */
void print_fluwfix(FILE*, char*, unsigned long, unsigned);

/* print text and unsigned long int in octal format to fix width */
void print_lowfix(char*, unsigned long, unsigned);
/* print text and unsigned long int in octal format to fix width to file */
void print_flowfix(FILE*, char*, unsigned long, unsigned);

/* print text and real number with fix line width */
void print_fwfix(char*, char*, double, unsigned);
/* print text and real number with fix line width to file */
void print_ffwfix(FILE*, char*, char*, double, unsigned);

/* print text and string to fix width */
void print_swfix(char*, char*, unsigned);
/* print text and string to fix width to file*/
void print_fswfix(FILE*, char*, char*, unsigned);

/* -------------------------------------------------------------------------- */

#endif
