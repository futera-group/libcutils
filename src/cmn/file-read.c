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
#include <stdlib.h>
#include "cmn/file.h"
#include "cmn/message.h"
#include "cmn/list.h"
#include "cmn/string.h"
#include "cmn/types.h"

/* -------------------------------------------------------------------------- */

/* add one value from column in file to the list
 
   q - the list
   t - data type
   w - pointer to value in string */
void file_read_column_add(struct list *q, short t, char *w) {
  char msg[] = "invalid format of input data in column";
  short unsigned d_su;
  long unsigned d_lu;
  unsigned d_u;
  double d_d;
  long double d_ld;
  short d_si;
  long d_li;
  int d_i;
  switch (t) {
    case TYPE_INT:
      if (sscanf(w,"%d",&d_i)!=1)
        msg_error(msg,1);
      list_add_end(q,&d_i,TYPE_INT);
      break;
    case TYPE_UINT:
      if (sscanf(w,"%u",&d_u)!=1)
        msg_error(msg,1);
      list_add_end(q,&d_u,TYPE_UINT);
      break;
    case TYPE_SINT:
      if (sscanf(w,"%hd",&d_si)!=1)
        msg_error(msg,1);
      list_add_end(q,&d_si,TYPE_SINT);
      break;
    case TYPE_USINT:
      if (sscanf(w,"%hu",&d_su)!=1)
        msg_error(msg,1);
      list_add_end(q,&d_su,TYPE_USINT);
      break;
    case TYPE_LINT:
      if (sscanf(w,"%ld",&d_li)!=1)
        msg_error(msg,1);
      list_add_end(q,&d_li,TYPE_LINT);
      break;
    case TYPE_ULINT: 
      if (sscanf(w,"%lu",&d_lu)!=1)
        msg_error(msg,1);
      list_add_end(q,&d_lu,TYPE_ULINT);
      break;
    case TYPE_DOUBLE:
      if (sscanf(w,"%lf",&d_d)!=1)
        msg_error(msg,1);
      list_add_end(q,&d_d,TYPE_DOUBLE);
      break;
    case TYPE_LDOUBLE:
      if (sscanf(w,"%Lf",&d_ld)!=1)
        msg_error(msg,1);
      list_add_end(q,&d_ld,TYPE_LDOUBLE);
      break;
    }
  }

/* read specified data column from file and return it in array
 
   f         - name of the file 
   col       - column in the file
   type      - data type
   row_first - first row for reading
   n_dat     - number of read data */
void *file_read_column(char *f, unsigned col, short type,
  long unsigned row_first, long unsigned *n_dat) {
  char *word,*line;
  long unsigned n_row;
  void *v = NULL;
  FILE *file;
  struct list *q;
  /* read data to fifo */
  file = file_open(f,"r");
  q = list_alloc();
  n_row = 0;
  line = str_read_line_new(file);
  while (line) {
    if (n_row>=row_first) {
      word = str_parse_word(line,col);
      file_read_column_add(q,type,word);
      }
    n_row++;
    line = str_free(line);
    line = str_read_line_new(file);
    }
  file_close(file);
  /* convert fifo to vector */
  v = list_conv_vec(q,type,n_dat);
  list_free(q,free);
  return(v);
  }

/* read specified data column interval from file and return it in array
 
   f         - name of the file 
   col       - column in the file
   type      - data type
   row_first - first row for reading
   row_last  - last row for reading
   n_dat     - number of read data */
void *file_read_column_nm(char *f, unsigned col, short type,
  long unsigned row_first, long unsigned row_last, long unsigned *n_dat) {
  char *line,*word;
  long unsigned n_row;
  void *v;
  FILE *file;
  struct list *q;
  /* read data to fifo */
  file = file_open(f,"r");
  q = list_alloc();
  n_row = 0;
  line = str_read_line_new(file);
  while (line) {
    if (n_row>=row_first) {
      word = str_parse_word(line,col);
      file_read_column_add(q,type,word);
      }
    line = str_free(line);
    n_row++;
    if (row_last && n_row>=row_last)
      break;
    line = str_read_line_new(file);
    }
  file_close(file);
  str_free(line);
  /* convert fifo to vector */
  v = list_conv_vec(q,type,n_dat);
  list_free(q,free);
  return(v);
  }

/* -------------------------------------------------------------------------- */
