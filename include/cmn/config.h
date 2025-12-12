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

#ifndef ZF_LIB_CMN_CONFIG_H
#define ZF_LIB_CMN_CONFIG_H

#include <stdio.h>
#include "cmn/list.h"
#include "cmn/queue.h"

/* -------------------------------------------------------------------------- */

/* config option */
struct cfg_item {
  char *key;                     /* keyword name */
  unsigned dep_level;            /* dependency level */
  short optional;                /* optionality indicator */
  short assigned;                /* assignment indicator */
  short data_type;               /* data type: int,float,double... */
  short data_rank;               /* data dimension: 1D,2D,3D... */
  void *data_ptr;                /* data value */
  void *data_dm1;                /* rank1 data dimensions */
  void *data_dm2;                /* rank2 data dimensions */
  struct list *item;             /* list of dependent items */
  void (*read_fce)(void*,void*,void*,FILE*); /* user-defined reading function */
  };

/* config file setting */
struct cfg {
  short verbose;                 /* verbose mode */
  short allow_free_lines;        /* free line allowance indicator */
  short allow_comment_lines;     /* comment line allowance indicator */
  short allow_comment_in_lines;  /* inline comments allowance indicator */
  short ignore_unknown_key;      /* unknown keywords ignoring */
  short ignore_unknown_frmt;     /* unknown line format ignoring */
  char comment_mark;             /* comment mark (default #) */
  char option_mark;              /* option mark (default $) */
  char dep_mark;                 /* dependency mark (default *) */
  char cont_mark;                /* continuation mark (default &) */
  struct list *item;             /* list of config-file items */
  };

/* -------------------------------------------------------------------------- */

/* allocate memory for new config file struct */
struct cfg *cfg_new(void);
/* allocate memory for new config item data struct */
struct cfg_item *cfg_item_new(void);

/* free memory allocated for config file struct */
void cfg_free(struct cfg*);
/* free memory of config-file item data saved in the list */
void cfg_item_free(struct cfg_item*);

/* check whether all mandatory options are assigned */
void cfg_check_assignment(struct cfg*);
/* search the item list and find element with specified keyword */
struct cfg_item* cfg_check_keyword(struct list*, char*, short);

/* define one option in config file and connect it with a variable */
void cfg_add_option(struct cfg*, char*, short, unsigned, short, void*,
  void*, void*, void(*fce)(void*,void*,void*,FILE*));
/* define dependent option in config file and connect it with a variable */
void cfg_add_option_dep(struct cfg*, char*, char*, short, unsigned, short,
  void*, void*, void*, void(*fce)(void*,void*,void*,FILE*));

/* assign data value from keyword line with specified type to variable */
void cfg_assign_rank0(char*, char*, short, void*);
/* assign data value from keyword line with specified type to variable */
void cfg_assign_rank1(struct queue*, char*, short, unsigned, void*);
/* assign data value from keyword line with specified type to variable */
void cfg_assign_rank2(struct queue**, char*, short, unsigned,
  unsigned*, void*);

/* allocate rank-1 data array in config item struct */
void cfg_alloc_rank1(void*, short, unsigned);
/* allocate rank-2 data array in config item struct */
void cfg_alloc_rank2(void*, short, unsigned);
/* allocate rank-3 data array in config item struct */
void cfg_alloc_rank3(void*, short, unsigned);
/* allocate rank-4 data array in config item struct */
void cfg_alloc_rank4(void*, short, unsigned);
/* allocate data array in config item struct */
void cfg_alloc_data(void*, short, int, unsigned);

/* return pointer to specific element in a rank-1 data array */
void *cfg_rank1_pointer(void*, short, unsigned, unsigned);
/* return pointer to specific element in a rank-2 data array */
void *cfg_rank2_pointer(void*, short, unsigned, unsigned);
/* return pointer to specific element in a rank-3 data array */
void *cfg_rank3_pointer(void*, short, unsigned, unsigned);
/* return pointer to specific element in a rank-4 data array */
void *cfg_rank4_pointer(void*, short, unsigned, unsigned);
/* return pointer to specific element in a data array */
void *cfg_data_pointer(void*, short, int, unsigned, unsigned);

/* define tree of sub-item templates in the new generated suboption */
void cfg_set_subsubitems(struct cfg*, struct cfg_item*, struct cfg_item*,
  unsigned, unsigned);
/* define subitems of the given options according to saved template */
void cfg_set_subitems(struct cfg*, struct cfg_item*);

/* read rank-1 data array from config file */
void cfg_read_data_rank1(struct cfg*, void*, short, unsigned, char*, FILE*);
/* read rank-2 data array from config file */
void cfg_read_data_rank2(struct cfg*, void*, short, void*, unsigned,
  char*, FILE*);
/* assign value from config file to keyword variable */
void cfg_read_keyword(struct cfg*, struct cfg_item*, char*, FILE*);
/* analyze keyword line from the config file */
void cfg_read_line(struct cfg*, char*, FILE*);
/* read config file */
void cfg_read(struct cfg*, char*);

/* -------------------------------------------------------------------------- */

#endif
