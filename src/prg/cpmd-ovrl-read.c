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
#include <cmn/file.h>
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include "prg/cpmd.h"

/* -------------------------------------------------------------------------- */

/* size of number types in CPMD binary files */
#define BIN_SIZE_INT   4  /* integer */
#define BIN_SIZE_REAL  8  /* double-precision real */

/* format of CPMD text files */
#define TXT_NUM_COLS   5  /* number of data columns */
#define TXT_VAL_WDTH  15  /* width reserved for data value */

/* -------------------------------------------------------------------------- */

/* read overlap matrix of AO orbitals from cpmd overlap file
 
   d    - main cpmd data struct
   file - name of the file */
void cpmd_ovrl_read_bin(struct cpmd_dat *d, char *file) {
  unsigned i,j;
  double x;
  int n;
  FILE *f;
  /* sanity check */
  if (sizeof(n)!=BIN_SIZE_INT)
    msg_error_f("inconsistent integer size for reading CPMD overlap matrix "
      "(%d/%d)",1,sizeof(n),BIN_SIZE_INT);
  if (sizeof(x)!=BIN_SIZE_REAL)
    msg_error_f("inconsistent real size for reading CPMD overlap matrix "
      "(%d/%d)",1,sizeof(x),BIN_SIZE_REAL);
  /* memory allocation */
  if (!d->ovrl)
    d->ovrl = mat_falloc(d->n_orbitals,d->n_orbitals);
  /* open file */
  f = file_open(file,"r");
  /* read data */
  if (fread(&n,BIN_SIZE_INT,1,f)!=1) /* record mark - begin */
    msg_error("CPMD overlap file is empty",1);
  if (n!=(BIN_SIZE_REAL*d->n_orbitals*d->n_orbitals))
    msg_error("inconsistent data record length in CPMD overlap file",1);
  for (i=0; i<d->n_orbitals; i++)
    for (j=0; j<d->n_orbitals; j++) {
      if (fread(&x,BIN_SIZE_REAL,1,f)!=1) /* overlap matrix element */
        msg_error_f("cannor read %d/%d element in"
          " CPMD overlap file",1,i+1,j+1);
      d->ovrl[i][j] = x;
      }
  if (fread(&n,BIN_SIZE_INT,1,f)!=1) /* record mark - end */
    msg_error("cannot read record end mark in CPMD overlap file",1);
  /* close file */
  file_close(f);
  }

/* read overlap matrix of AO orbitals from text-formatted cpmd overlap file
 
   d    - main cpmd data struct
   file - name of the file */
void cpmd_ovrl_read_txt(struct cpmd_dat *d, char *file) {
  char *line;
  unsigned i,j,k,n,n_blocks,n_rows,n_cols;
  FILE *f;
  /* memory allocation */
  if (!d->ovrl)
    d->ovrl = mat_falloc(d->n_orbitals,d->n_orbitals);
  /* open file */
  f = file_open(file,"r");
  /* number of data blocks */
  n_blocks = d->n_orbitals/TXT_NUM_COLS;
  if (d->n_orbitals%TXT_NUM_COLS)
    n_blocks++;
  n_rows = d->n_orbitals;
  /* data blocks */
  for (i=0; i<n_blocks; i++) {
    /* header */
    line = str_read_line_new(f);
    if (!line)
      msg_error_f("cannot read head of data block #%d in"
        " CPMD overlap file",1,i+1);
    line = str_free(line);
    /* data rows */
    for (j=0; j<n_rows; j++) {
      line = str_read_line_new(f);
      if (!line)
        msg_error_f("cannot read data row #%d of block #%d in"
          " CPMD overlap file",1,j+1,i+1);
      /* number of columns */
      n_cols = (j<TXT_NUM_COLS ? j+1 : TXT_NUM_COLS);
      /* offset */
      n = 0;
      while (n<str_length(line) && line[n]==' ') n++;
      while (n<str_length(line) && line[n]!=' ') n++;
      /* coefficients */
      for (k=0; k<n_cols; k++)
        if (sscanf(line+n+k*TXT_VAL_WDTH,"%lf",
          &(d->ovrl[d->n_orbitals-n_rows+j][i*TXT_NUM_COLS+k]))!=1)
          msg_error_f("cannot read %d/%d element in CPMD overlap file",
            1,d->n_orbitals-n_rows+j+1,i*TXT_NUM_COLS+k+1);
      line = str_free(line);
      }
    /* number of rows */
    n_rows -= TXT_NUM_COLS;
    }
  /* close file */
  file_close(f);
  }

/* -------------------------------------------------------------------------- */
