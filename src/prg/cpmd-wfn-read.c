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

#include <cmn/file.h>
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

/* read projection of states from cpmd wavefunction file
 
   d    - main cpmd data struct
   file - name of the file */
void cpmd_wfn_read_bin(struct cpmd_dat *d, char *file) {
  unsigned i,j;
  double x;
  int n;
  FILE *f;
  /* sanity check */
  if (sizeof(n)!=BIN_SIZE_INT)
    msg_error_f("inconsistent integer size for reading CPMD wavefunction "
      "(%d/%d)",1,sizeof(n),BIN_SIZE_INT);
  if (sizeof(x)!=BIN_SIZE_REAL)
    msg_error_f("inconsistent real size for reading CPMD wavefunction "
      "(%d/%d)",1,sizeof(x),BIN_SIZE_REAL);
  /* open file */
  f = file_open(file,"r");
  /* read header */
  if (fread(&n,BIN_SIZE_INT,1,f)!=1) /* record mark - begin */
    msg_error("CPMD wavefunction file is empty",1);
  if (fread(&n,BIN_SIZE_INT,1,f)!=1) /* number of AOs */
    msg_error("cannor read number of AOs in CPMD wavefunction file",1);
  if (n!=d->n_orbitals)
    msg_error_f("inconsistent number of AOs in CPMD wavefunction file"
      " (%d/%d)",1,n,d->n_orbitals);
  if (fread(&n,BIN_SIZE_INT,1,f)!=1) /* number of types */
    msg_error("cannor read number of atom types in CPMD wavefunction file",1);
  if (n!=d->n_types)
    msg_error_f("inconsistent number of atom types in CPMD wavefunction file"
      " (%d/%d)",1,n,d->n_types);
  for (i=0; i<d->n_types; i++) {
    if (fread(&x,BIN_SIZE_REAL,1,f)!=1) /* valence charge */
      msg_error_f("cannor read valence charge of type #%d in"
        " CPMD wavefunction file",1,i+1);
    if (fread(&n,BIN_SIZE_INT,1,f)!=1)  /* number of atoms */
      msg_error_f("cannor read number of type #%d atoms in"
        " CPMD wavefunction file",1,i+1);
    if (fread(&n,BIN_SIZE_INT,1,f)!=1)  /* number of AOs */
      msg_error_f("cannor read number of type #%d AOs in"
        " CPMD wavefunction file",1,i+1);
    }
  if (fread(&n,BIN_SIZE_INT,1,f)!=1) /* record mark - end */
    msg_error("cannot read header end mark in CPMD wavefunction file",1);
  /* read wavefunction */
  if (fread(&n,BIN_SIZE_INT,1,f)!=1) /* record mark - begin */
    msg_error("cannot read wavefunction record in CPMD wavefunction file",1);
  if (n!=(BIN_SIZE_REAL*d->n_states*d->n_orbitals))
    msg_error("inconsistent wavefunction record length in "
      "CPMD wavefunction file",1);
  for (i=0; i<d->n_states; i++)
    for (j=0; j<d->n_orbitals; j++) {
      if (fread(&x,BIN_SIZE_REAL,1,f)!=1) /* wfn coeff */
        msg_error_f("cannor read coefficient %d/%d in"
          " CPMD wavefunction file",1,i+1,j+1);
        d->state_a[i].ao_pj[j] = x;
        }
  if (fread(&n,BIN_SIZE_INT,1,f)!=1) /* record mark - end */
    msg_error("cannot read wavefunction end mark in CPMD wavefunction file",1);
  /* close file */
  file_close(f);
  }

/* read projection of states from cpmd text-formatted wavefunction file
 
   d    - main cpmd data struct
   file - name of the file */
void cpmd_wfn_read_txt(struct cpmd_dat *d, char *file) {
  char *line;
  unsigned i,j,k,n,n_blocks,n_cols;
  FILE *f;
  /* open file */
  f = file_open(file,"r");
  /* number of data blocks */
  n_blocks = d->n_states/TXT_NUM_COLS;
  if (d->n_states%TXT_NUM_COLS)
    n_blocks++;
  /* data blocks */
  for (i=0; i<n_blocks; i++) {
    /* number of columns */
    n_cols = ((i+1)<n_blocks || !(d->n_states/TXT_NUM_COLS) ?
      TXT_NUM_COLS : d->n_states%TXT_NUM_COLS);
    /* header */
    line = str_read_line_new(f);
    if (!line)
      msg_error_f("cannot read head of data block #%d in"
        " CPMD wavefunction file",1,i+1);
    line = str_free(line);
    /* data rows */
    for (j=0; j<d->n_orbitals; j++) {
      line = str_read_line_new(f);
      if (!line)
        msg_error_f("cannot read data row #%d of block #%d in"
          " CPMD wavefunction file",1,j+1,i+1);
      /* offset */
      n = 0;
      while (n<str_length(line) && line[n]==' ') n++;
      while (n<str_length(line) && line[n]!=' ') n++;
      /* coefficients */
      for (k=0; k<n_cols; k++)
        if (sscanf(line+n+k*TXT_VAL_WDTH,"%lf",
          &(d->state_a[i*TXT_NUM_COLS+k].ao_pj[j]))!=1)
          msg_error_f("cannot read AO coefficient #%d of state #%d in"
            " CPMD wavefunction file",1,j+1,i*TXT_NUM_COLS+k);
      line = str_free(line);
      }
    }
  /* close file */
  file_close(f);
  }

/* -------------------------------------------------------------------------- */
