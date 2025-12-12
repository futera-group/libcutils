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

#ifndef ZF_LIB_CMN_MESSAGE_H
#define ZF_LIB_CMN_MESSAGE_H

/* -------------------------------------------------------------------------- */

/* set pointer to error function */
void msg_error_fce_set(void(*e_fce)(void*,void*,void*), void*, void*,void*);

/* print error message */
void msg_error(char*, short);
/* print formatted error message */
void msg_error_f(char*, short, ...);

/* print warning */
void msg_warn(char*);
/* print warning with user specified format */
void msg_warn_f(char*, ...);

/* print debug message */
void msg_debug(char*, short);
/* print debug message and flusth the output stream */
void msg_Debug(char*);
/* print function debug message and flush the output stream */
void msg_Debug_fce(char*, char*);

/* -------------------------------------------------------------------------- */

#endif
