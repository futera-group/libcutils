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
#include "cmn/config.h"
#include "cmn/file.h"
#include "cmn/message.h"
#include "cmn/queue.h"
#include "cmn/string.h"
#include "cmn/types.h"
#include "cmn/vector.h"

/* -------------------------------------------------------------------------- */

/* read rank-1 data array from config file

   c    - pointer to config data struct
   d    - pointer to the data storage
   type - type of the data
   nd   - number of expected values
   key  - config-file keyword
   f    - open file stream */
void cfg_read_data_rank1(struct cfg *c, void *d, short type,
  unsigned nd, char *key, FILE *f) {
  char *line;
  struct queue *q;
  /* read data from file */
  q = queue_alloc();
  if (nd) {
    for (;;) {
      line = str_read_line_new(f);
      if (!line)
        msg_error_f("unexpected end of file while reading \"%s\" data",1,key);
      str_trim(line);
      if (c->allow_comment_in_lines)
        str_rtrim_mark(line,c->comment_mark);
      if (str_char_last_nb(line,NULL)!=c->cont_mark) {
        queue_sadd(q,line);
        break;
        }
      else {
        str_rtrim_mark(line,c->cont_mark);
        str_rtrim_mark(line,',');
        queue_sadd(q,line);
        }
      str_free(line);
      }
    }
  /* allocate data storage */
  cfg_alloc_rank1(d,type,nd);
  /* data assignment */
  cfg_assign_rank1(q,key,type,nd,d);
  /* clean memory */
  queue_free(q);
  }

/* read rank-2 data array from config file

   c    - pointer to config data struct
   d    - pointer to the data storage
   type - type of the data
   dm   - storage for rank-1 array dimensions
   nd   - number of expected values
   key  - config-file keyword
   f    - open file stream */
void cfg_read_data_rank2(struct cfg *c, void *d, short type,
  void *dm, unsigned nd, char *key, FILE *f) {
  unsigned i;
  char *line;
  struct queue **q;
  /* read data from file */
  q = vec_talloc(sizeof(struct queue*),nd);
  for (i=0; i<nd; i++) {
    q[i] = queue_alloc();
    for (;;) {
      line = str_read_line_new(f);
      if (!line)
        msg_error_f("unexpected end of file while reading \"%s\" data",1,key);
      str_trim(line);
      if (c->allow_comment_in_lines)
        str_rtrim_mark(line,c->comment_mark);
      queue_sadd(q[i],line);
      if (str_char_last_nb(line,NULL)!=c->cont_mark)
        break;
      str_free(line);
      }
    }
  /* allocate data storage */
  cfg_alloc_rank1(dm,TYPE_UINT,nd);
  cfg_alloc_rank2(d,type,nd);
  /* data assignment */
  cfg_assign_rank2(q,key,type,nd,*(unsigned**)dm,d);
  /* clean memory */
  for (i=0; i<nd; i++)
    queue_free(q[i]);
  vec_tfree(q);
  }

/* assign value from config file to keyword variable

   d - pointer to config data struct
   c - pointer to config item struct
   s - line from config file
   f - open file stream */
void cfg_read_keyword(struct cfg *d, struct cfg_item *c, char *s, FILE *f) {
  /* user-defined reading */
  if (c->read_fce)
    cfg_assign_rank0(s,c->key,c->data_type,c->data_dm1);
  /* standard reading */
  else {
    switch (c->data_rank) {
      case 0: /* scalar value */
        cfg_assign_rank0(s,c->key,c->data_type,c->data_ptr);
        break;
      case 1: /* 1D array */
        cfg_assign_rank0(s,c->key,TYPE_UINT,c->data_dm1);
        cfg_read_data_rank1(d,c->data_ptr,c->data_type,
          *((unsigned*)c->data_dm1),c->key,f);
        break;
      case 2: /* 2D array */
        cfg_assign_rank0(s,c->key,TYPE_UINT,c->data_dm2);
        cfg_read_data_rank2(d,c->data_ptr,c->data_type,c->data_dm1,
          *((unsigned*)c->data_dm2),c->key,f);
        break;
      default:
        msg_error_f("unsupported data rank (%d) of \"%s\" keyword",1,
          c->data_rank,c->key);
      }
    }
  /* dependent options */
  cfg_set_subitems(d,c);
  /* user-defined data */
  if (c->read_fce)
    c->read_fce(c->data_ptr,c->data_dm1,c->data_dm2,f);
  }

/* analyze keyword line from the config file

   c - pointer to config file setting struct
   s - line from config file
   f - pointer to open config file */
void cfg_read_line(struct cfg *c, char *s, FILE *f) {
  char key[1024],frmt[10];
  struct cfg_item *p;
  /* unknown data line */
  if (s[0]!=c->option_mark) {
    if (!c->ignore_unknown_frmt)
      msg_error_f("unexpected data: %s",1,s);
    else if (c->verbose)
      printf("unexpected format, ignored...\n");
    }
  /* read option */
  sprintf(frmt,"%c%%s",c->option_mark);
  sscanf(s,frmt,key);
  /* unknown option */
  p = cfg_check_keyword(c->item,key,0);
  if (!p) {
    if (!c->ignore_unknown_key)
      msg_error_f("unknown keyword \"%s\"",1,key);
    else if (c->verbose)
      printf("unknown keyword \"%s\", ignored...\n",key);
    }
  /* assign value */
  else {
    cfg_read_keyword(c,p,s,f);
    p->assigned=1;
    }
  }

/* read config file 

   c - pointer to config file setting
   f - config file name */
void cfg_read(struct cfg *c, char *f) {
  char *line;
  FILE *file;
  /* open the file */
  file = file_open(f,"r");
  line = str_read_line_new(file);
  while (line) {
    /* trim whites */
    str_trim(line);
    /* trim comments */
    if (c->allow_comment_in_lines)
      str_rtrim_mark(line,c->comment_mark);
    /* free line */
    if (str_length(line)<2) {
      if (!c->allow_free_lines)
        msg_error("unexpected free line",1);
      }
    /* commented lines */
    else if (c->allow_comment_lines && line[0]==c->comment_mark) {
      }
    /* read setting */
    else
      cfg_read_line(c,line,file);
    /* next line */
    line = str_free(line);
    line = str_read_line_new(file);
    }
  /* check whether all mandatory options are set */
  cfg_check_assignment(c);
  /* close the file */
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
