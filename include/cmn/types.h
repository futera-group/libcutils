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

#ifndef ZF_LIB_CMN_TYPES_H
#define ZF_LIB_CMN_TYPES_H

/* -------------------------------------------------------------------------- */

/* symbolic constants */
#define TYPE_VOID    0  /* void */
#define TYPE_INT     1  /* integer */
#define TYPE_UINT    2  /* unsigned integer */
#define TYPE_SINT    3  /* short integer */
#define TYPE_USINT   4  /* unsigned short integer */
#define TYPE_LINT    5  /* long integer */
#define TYPE_ULINT   6  /* unsigned long integer */
#define TYPE_DOUBLE  7  /* double precision real */
#define TYPE_LDOUBLE 8  /* long double precision real */
#define TYPE_STRING  9  /* string of chars */

/* -------------------------------------------------------------------------- */

/* return string name of the specified data type */
char* type_name(short);
/* return string flag of the specified data type */
char* type_flag(short);

/* -------------------------------------------------------------------------- */

/* allocate memory for one specified data type */
void type_alloc(void*, short);
/* allocate memory for vector of specified data types */
void type_alloc_v(void*, short, unsigned, ...);
/* free memory allocated by type_alloc functions */
void* type_free(void*);

/* allocate and copy value of of a given data type */
void *type_copy(void*, short);

/* read value of specified type from string */
short type_read(char*, void*, short);

/* put one value into array of specified type */
void type_array_put(void*, unsigned, void*, short);

/* -------------------------------------------------------------------------- */

#endif
