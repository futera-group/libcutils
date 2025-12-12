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

#ifndef ZF_LIB_CMN_STRING_H
#define ZF_LIB_CMN_STRING_H

#include <stdio.h>

/* -------------------------------------------------------------------------- */

#define STR_BUFFER_LENGTH 4098

/* -------------------------------------------------------------------------- */

/* basic operations */

/* allocate memory for a string */
char* str_new(unsigned);
/* free memory allocated for string */
char* str_free(char*);

/* return lenght of a string */
unsigned str_length(char*);

/* copy a string */
void str_copy(char*, char*);
/* creates new copy of the character string */
char* str_copy_new(char*);

/* merge two strings */
void str_merge(char*, char*, char*);
/* merge two string to a new one */
char* str_merge_new(char*, char*);

/* append string to dynamically allocated one */
void str_append(char**, char*);

/* formatted print to new dynamically-allocated string */
char* str_sprintf(char*, ...);

/* function operating with file names or pathes */

/* cut name of file from the path without extension */
int str_file_base(char*, char*);
/* return base nae of a file */
char* str_file_base_new(char*);

/* cut name of file from the path */
int str_file_name(char*, char*);
/* cut name of file from the path and return it */
char* str_file_name_new(char*);

/* cut file name extension from the path */
int str_file_ext(char*, char*);
/* return extension of a file */
char* str_file_ext_new(char*);

/* cut name of directory from the path */
int str_file_dir(char*, char*);
/* cut name of directory from the path and return it */
char* str_file_dir_new(char*);

/* create full path from directory and file name */
int str_file_path(char*, char*, char*);
/* return full path to a file */
char *str_file_path_new(char*, char*);

/* function changing letters in string (uppercase <-> lowercase) */

/* convert string to lowcase chars, the string is changed */
char* str_lowcase(char*);
/* convert string to lowcase chars, the string is not changed */
void str_lowcase_copy(char*, char*);
/* create new lowercase version of input string */
char* str_lowcase_new(char*);

/* convert string to lowcase chars instead of the first which is capital */
char* str_Lowcase(char*);
/* convert string to lowcase chars instead of the first which is capital */
void str_Lowcase_copy(char*, char*);
/* create new lowercase version of input string with first capital */
char* str_Lowcase_new(char*);

/* convert string to upcase chars, the string is changed */
char* str_upcase(char*);
/* convert string to upcase chars, the string is not changed */
void str_upcase_copy(char*, char*);
/* create new upper-case version of input string */
char* str_upcase_new(char*);

/* functions for trimming */

/* trim spaces from both sides of the string, the string is changed */
void str_trim(char*);
/* trim spaces from both sides of the string, the string is not changed */
void str_trim_copy(char*, char*);
/* trim string from both sides up to specified character */
void str_trim_mark(char*, char);
/* delete specific characters from the string */
void str_trim_char(char*, char*, char*);

/* trim spaces from right side of the string */
void str_rtrim(char*);
/* trim string from right side up to specified character */
void str_rtrim_mark(char*, char);

/* trim spaces from left side of the string */
void str_ltrim(char*);
/* trim string from left side up to specified character */
void str_ltrim_mark(char*, char);

/* functions for parsing and splitting */

/* split string to words according to given deliminer */
char ** str_split(char*, char, unsigned*);
/* split string to words according to given string deliminer */
char** str_split_s(char*, char*, unsigned*);
/* split string to words according to given word deliminer */
char** str_split_w(char*, char*, unsigned*);

/* convert string to array of numbers with specified type */
void str_parse_array(char*, short, void*, unsigned*);
/* convert string to array of int numbers */
int *str_parse_iarray(char*, unsigned*);
/* convert string to array of u-int numbers */
unsigned *str_parse_uarray(char*, unsigned*);
/* convert string to array of s-int numbers */
short *str_parse_siarray(char*, unsigned*);
/* convert string to array of su-int numbers */
short unsigned *str_parse_suarray(char*, unsigned*);
/* convert string to array of l-int numbers */
long *str_parse_liarray(char*, unsigned*);
/* convert string to array of lu-int numbers */
long unsigned *str_parse_luarray(char*, unsigned*);
/* convert string to array of real numbers */
double *str_parse_farray(char*, unsigned*);
long double *str_parse_lfarray(char*, unsigned*);

/* parse string containing list of numbers separated by commas and dashes */
void str_parse_num(char*, short, void*, unsigned*);
/* parse general string array with u-int numbers */
unsigned *str_parse_unum(char*, unsigned*);

/* read n-th word from the string and return pointer to it */
char *str_parse_word(char*, unsigned);

/* functions dealing with numbers / characters */

/* return lenght of string representing int */
unsigned str_inum_len(int);
/* return length of string representing unsigned int */
unsigned str_unum_len(unsigned);
/* return lenght of string representing int */
unsigned str_sinum_len(short);
/* return lenght of string representing unsigned short */
unsigned str_sunum_len(unsigned short);
/* return lenght of string representing long int */
unsigned str_linum_len(long);
/* return length of string representing unsigned long */
unsigned str_lunum_len(unsigned long);

/* convert vector of numbers to string representation (u-int) */
void str_unum_range(char*, unsigned*, unsigned);
/* convert vector of c-indices to string representation (u-int) */
void str_unum_Range(char*, unsigned*, unsigned);

/* convert fortran string number to double precision real */
double str_fnum_d2e(char*);

/* convert unsigned integer number to string representation */
char* str_unum_i2s_new(unsigned);

/* count number of specific characters in given string */
unsigned str_char_count(char*, char);
/* get the first blank character in the string */
char str_char_first_b(char*, unsigned*);
/* get the first non-blank character in the string  */
char str_char_first_nb(char*, unsigned*);
/* get last non-blank character in the string */
char str_char_last_nb(char*, unsigned*);

/* function treating with substrings */

/* copy substring interval */
void str_sub_copy(char*, char*, unsigned, unsigned);
/* copy substring after specified deliminer */
void str_sub_copy_adel(char*, char*, char);
/* copy substring before specified deliminer */
void str_sub_copy_bdel(char*, char*, char);
/* find substring at the begin of other string */
int str_sub_bfind(char*, char*);

/* functions for reading strings */

/* read line form input stream (unspecified lenght) */
char* str_read_line_new(FILE *);

/* read line from standard input and return first word */
void str_read_word(char*);
/* read line from standard input and return new string with first word */
char* str_read_word_new(void);
/* read one word from the given string */
unsigned str_read_word_s(char*, char*, char, unsigned);

/* read yes or no answer to given question */
short str_read_yesno(char*);
/* read yes or no answer to given question with pre-set default */
short str_read_yesno_def(char*, short);

/* ask for and integer number */
int str_read_inum(char*);
/* ask for and read unsigned integer number */
unsigned str_read_unum(char*);
/* ask for and read unsigned integer number, offer default */
unsigned str_read_unum_def(char*, unsigned);
/* ask for and read unsigned integer number from given interval */
unsigned str_read_unum_int(char*, unsigned, unsigned);
/* ask for and read double precision real number */
double str_read_fnum(char*);

/* string comparisons */

/* compare two strings */
int str_compare(char*, char*);
/* case-insensitive comparison of two strings */
int str_compare_nc(char*, char*);

/* string search */

/* find string in a file */
int str_ffind(FILE*, char*, char*);
/* find word in a file, stop at first occurence and return the line */
char* str_ffind_new(FILE*, char*);
/* find word in a file, at the beginning of line, stop at first occurrence */
int str_ffind_b(FILE*, char*, char*);
/* find word in a file, at the beggining of line, return the line */
char* str_ffind_b_new(FILE*, char*);

/* substitute all specified characters in the string */
void str_subst(char*, char, char);
/* substitute first specified characters in the string */
void str_subst_one(char*, char, char);
/* substitute first specified characters in the string by u-int */
void str_subst_one_ui(char**, char, unsigned);

/* -------------------------------------------------------------------------- */

#endif
