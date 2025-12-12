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
#include "cmn/file.h"
#include "cmn/message.h"

/* buffer size for reading */
#define FILE_BIN_CMAX 10000

/* -------------------------------------------------------------------------- */

/* copy n bytes from one binary file to another

   file1 - the file from which data are read
   file2 - the file to which data are writen
   nbyte - number of bytes */
void file_bin_copy(FILE *file1, FILE *file2, unsigned long nbyte) {
  unsigned long nrec;
  static char cstore[FILE_BIN_CMAX];
  while (nbyte>0) {
    nrec = (nbyte>FILE_BIN_CMAX ? FILE_BIN_CMAX : nbyte);
    if (fread(&cstore,nrec,1,file1)<1)
      msg_error("reading error during copying binary file",1);
    if (fwrite(&cstore,nrec,1,file2)<1)
      msg_error("writing error during copying binary file",1);
    nbyte-=nrec;
    }
  }

/* copy n bytes from one binary file to another and save them in array

   file1 - the file from which data are read
   file2 - the file to which data are writen
   nitem - number of items in the array
   nsize - size of one item in the array in bytes
   store - the array for saved data */
void file_bin_copy_save(FILE *file1, FILE *file2, unsigned long nitem,
  unsigned nsize, void *store) {
  if (fread(store,nsize,nitem,file1)<nitem)
    msg_error("reading error during copying binary file",1);
  if (fwrite(store,nsize,nitem,file2)<nitem)
    msg_error("writing error during copying binary file",1);
  }

/* copy rest of the binary file to another file 
  
   file1 - the file from which data are read
   file2 - the file to which data are writen */
void file_bin_copy_all_f(FILE *file1, FILE *file2) {
  char c;
  while (fread(&c,1,1,file1)>0)
    if (fwrite(&c,1,1,file2)<1)
      msg_error("writing error during copying binary file",1);
  }

/* copy rest of the binary file to another file 
  
   file1 - the file from which data are read
   file2 - the file to which data are writen */
void file_bin_copy_all(char *file1, char *file2) {
  FILE *f1,*f2;
  f1 = file_open(file1,"r");
  f2 = file_open(file2,"w");
  file_bin_copy_all_f(f1,f2);
  file_close(f1);
  file_close(f2);
  }

/* -------------------------------------------------------------------------- */

/* Read string of characters from Fortran binary file

   d - storage for the data
   n - number of characters
   f - open file stream */
void file_bin_frt_read_str(char *d, unsigned n, FILE *f) {
  int rec_l;
  char c;
  /* sanity check */
  if (sizeof(rec_l)!=FILE_BIN_FRT_SIZE_INT)
    msg_error_f("incompatible integer size (%d) for fortran data reading",
      1,sizeof(rec_l));
  if (sizeof(c)!=FILE_BIN_FRT_SIZE_CHAR)
    msg_error_f("incompatible character size (%d) for fortran data reading",
      1,sizeof(c));
  /* header */
  if (fread(&rec_l,FILE_BIN_FRT_SIZE_INT,1,f)<1)
    msg_error("cannot read fortran record length",1);
  /* values */
  if (fread(d,FILE_BIN_FRT_SIZE_CHAR,n,f)<1)
    msg_error("cannot read fortran string record",1);
  /* footer */
  }

/* Read array of integers from Fortran binary file

   d - storage for the data
   n - number of values
   f - open file stream */
void file_bin_frt_read_int(int *d, unsigned n, FILE *f) {
  int rec_l;
  /* sanity check */
  if (sizeof(rec_l)!=FILE_BIN_FRT_SIZE_INT)
    msg_error_f("incompatible integer size (%d) for fortran data reading",
      1,sizeof(rec_l));
  /* header */
  if (fread(&rec_l,FILE_BIN_FRT_SIZE_INT,1,f)<1)
    msg_error("cannot read fortran record length",1);
  /* values */
  if ((n*FILE_BIN_FRT_SIZE_INT)!=rec_l)
    msg_error_f("incompatible fortran integer-data record length (%d/%d)",
      1,rec_l,n*FILE_BIN_FRT_SIZE_INT);
  if (fread(d,FILE_BIN_FRT_SIZE_INT,n,f)<1)
    msg_error("cannot read fortran integer data record",1);
  /* footer */
  if (fread(&rec_l,FILE_BIN_FRT_SIZE_INT,1,f)<1)
    msg_error("cannot read fortran record ending",1);
 }

/* Read array of real numbers from Fortran binary file

   d - storage for the data
   n - number of values
   f - open file stream */
void file_bin_frt_read_real(double *d, unsigned n, FILE *f) {
  int rec_l;
  double x;
  /* sanity check */
  if (sizeof(rec_l)!=FILE_BIN_FRT_SIZE_INT)
    msg_error_f("incompatible integer size (%d) for fortran data reading",
      1,sizeof(rec_l));
  if (sizeof(x)!=FILE_BIN_FRT_SIZE_REAL)
    msg_error_f("incompatible real-number size (%d) for fortran data reading",
      1,sizeof(x));
  /* header */
  if (fread(&rec_l,FILE_BIN_FRT_SIZE_INT,1,f)<1)
    msg_error("cannot read fortran record length",1);
  /* values */
  if ((n*FILE_BIN_FRT_SIZE_REAL)!=rec_l)
    msg_error_f("incompatible fortran real-data record length (%d/%d)",
      1,rec_l,n*FILE_BIN_FRT_SIZE_REAL);
  if (fread(d,FILE_BIN_FRT_SIZE_REAL,n,f)<1)
    msg_error("cannot read fortran real data record",1);
  /* footer */
  if (fread(&rec_l,FILE_BIN_FRT_SIZE_INT,1,f)<1)
    msg_error("cannot read fortran record ending",1);
  }

/* Read array of single-precsion real numbers from Fortran binary file

   d - storage for the data
   n - number of values
   f - open file stream */
void file_bin_frt_read_real_s(float *d, unsigned n, FILE *f) {
  int rec_l;
  float x;
  /* sanity check */
  if (sizeof(rec_l)!=FILE_BIN_FRT_SIZE_INT)
    msg_error_f("incompatible integer size (%d) for fortran data reading",
      1,sizeof(rec_l));
  if (sizeof(x)!=FILE_BIN_FRT_SIZE_REAL_S)
    msg_error_f("incompatible s-real-number size (%d) for fortran data reading",
      1,sizeof(x));
  /* header */
  if (fread(&rec_l,FILE_BIN_FRT_SIZE_INT,1,f)<1)
    msg_error("cannot read fortran record length",1);
  /* values */
  if ((n*FILE_BIN_FRT_SIZE_REAL_S)!=rec_l)
    msg_error_f("incompatible fortran s-real-data record length (%d/%d)",
      1,rec_l,n*FILE_BIN_FRT_SIZE_REAL_S);
  if (fread(d,FILE_BIN_FRT_SIZE_REAL_S,n,f)<1)
    msg_error("cannot read fortran real data record",1);
  /* footer */
  if (fread(&rec_l,FILE_BIN_FRT_SIZE_INT,1,f)<1)
    msg_error("cannot read fortran record ending",1);
  }

/* -------------------------------------------------------------------------- */
