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

#ifndef ZF_LIB_CMN_DIR_H
#define ZF_LIB_CMN_DIR_H

#include <dirent.h>
#include <sys/types.h>

/* -------------------------------------------------------------------------- */

/* detect if given directory exists */
int dir_exist(char*);

/* get the current directory path */
void dir_current(char*);
/* return new string with current working directory path */
char* dir_current_new(void);

/* create directory with given mask */
int dir_make(char*, mode_t, short);

/* change current working directory */
int dir_change(char*, short);

/* open directory and return its pointer */
DIR *dir_open(char*, short);
/* open the current directory for reading */
DIR *dir_open_current(void);

/* close opened directory */
void dir_close(DIR*);

/* -------------------------------------------------------------------------- */

#endif
