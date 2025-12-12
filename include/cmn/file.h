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

#ifndef ZF_LIB_CMN_FILE_H
#define ZF_LIB_CMN_FILE_H

#include <sys/stat.h>
#include <stdio.h>

/* fortran files */
#define FILE_BIN_FRT_SIZE_CHAR    1  /* size of character */
#define FILE_BIN_FRT_SIZE_INT     4  /* size of integer */
#define FILE_BIN_FRT_SIZE_REAL    8  /* size of double */
#define FILE_BIN_FRT_SIZE_REAL_S  4  /* size of float */

/* -------------------------------------------------------------------------- */

/* file data struct */
struct file_dat {
  char *name;     /* file name */
  char *path;     /* path to the file */
  };

/* -------------------------------------------------------------------------- */

/* create new file info data struct from given path */
struct file_dat *file_dat_new(char*);
/* free memory allocated for the file data struct */
void file_dat_free(struct file_dat*);

/* open file and return its pointer */
FILE *file_open(char*, char*);
/* open temporary file */
FILE *file_open_tmp(void);
/* associate open file stream with different file */
FILE *file_reopen(char*, char*, FILE*);
/* close opened file */
void file_close(FILE*);

/* get information about file */
int file_stat(char*, struct stat*, short);
/* detect if given file exists */
int file_exist(char*);
/* return size of specified file */
long unsigned file_size(char*);

/* create symbolic link to file saved in the filesystem */
void file_link(char*, char*);
/* copy one file to another */
void file_copy(char*, char*);
/* move or rename file saved in the specified path */
void file_move(char*, char*);
/* merge two files to new one */
void file_merge(char*, char*, char*);
/* remove a file */
void file_remove(char*);

/* find path to given executable file */
void file_exec_path(char*);
/* return string path to given executable file */
char* file_exec_path_new(char*);

/* read specified data column from file and return it in array */
void *file_read_column(char*, unsigned, short, long unsigned, long unsigned*);
/* read specified data column interval from file and return it in array */
void *file_read_column_nm(char*, unsigned, short, long unsigned, long unsigned,
  long unsigned*);

/* search and create list of all file of specified type */
struct list *file_search(char*, short, char*);

/* binary file function */

/* copy n bytes from one binary file to another */
void file_bin_copy(FILE*, FILE*, unsigned long);
/* copy n bytes from one binary file to another and save them in array*/
void file_bin_copy_save(FILE*, FILE*, unsigned long, unsigned, void*);
/* copy rest of the binary file to another file */
void file_bin_copy_all_f(FILE*, FILE*);
/* copy rest of the binary file to another file */
void file_bin_copy_all(char*, char*);

/* fortran binary files */

/* Read string of characters from Fortran binary file */
void file_bin_frt_read_str(char*, unsigned, FILE*);
/* read array of integers from Fortran binary file */
void file_bin_frt_read_int(int*, unsigned, FILE*);
/* read array of real numbers from Fortran binary file */
void file_bin_frt_read_real(double*, unsigned, FILE*);
/* Read array of single-precsion real numbers from Fortran binary file */
void file_bin_frt_read_real_s(float*, unsigned, FILE*);

/* -------------------------------------------------------------------------- */

#endif
