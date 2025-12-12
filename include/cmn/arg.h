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

#ifndef ZF_LIB_CMN_ARG_H
#define ZF_LIB_CMN_ARG_H

#include "cmn/list.h"

/* -------------------------------------------------------------------------- */

/* structs and types */
struct arg_dat {
  char *name;               /* one letter representing the option */
  char *desc;               /* short description of the option */
  short type;               /* type of the assigned argument */
  short swtch;              /* switch not option, no argument needed */
  short optional;           /* option is optional not mandatory */
  short simple;             /* simple argument without specifier */
  short set;                /* specify if the argument was already read */
  struct arg_dat *bind;     /* specify set of bound options */
  struct arg_dat *excl;     /* specify set of excluded options */
  void *val;                /* place for saving argument value */
  };

struct arg {
  short help;
  void (*help_fce)(void);
  struct list *data;
  };

/* -------------------------------------------------------------------------- */

/* allocate memory for the command-line argument structure */
struct arg* arg_new(void);
/* allocate memory for argument data structures */
struct arg_dat* arg_dat_new(void);

/* free memory allocated for a command-line argument structure */
struct arg* arg_free(struct arg*);
/* free memory allocated for a command-line argument data structure */
struct arg_dat* arg_dat_free(struct arg_dat*);

/* define switch */
struct arg_dat* arg_add_switch(struct arg*, char*, short, short*);
/* define option */
struct arg_dat* arg_add_option(struct arg*, char*, char*, short, short, void*);
/* define simple argument */
struct arg_dat* arg_add_item(struct arg*, char*, short, short, void*);
/* assign help function */
void arg_add_help(struct arg*, void (*fce)(void));

/* set option binding */
void arg_bind(struct arg*, struct arg_dat*, struct arg_dat*);
/* set option exclusion */
void arg_exclude(struct arg*, struct arg_dat*, struct arg_dat*);

/* read and analyze command line arguments */
void arg_read(struct arg*, unsigned, char**);

/* print usage of the program */
void arg_print_usage(struct arg*, char*, char*, short);

/* convert string to the correct type and assign it with the variable */
int arg_var_set(char*, void*, short);

/* check presence and compatibility */
void arg_check(struct arg*, char*);

/* -------------------------------------------------------------------------- */

#endif
