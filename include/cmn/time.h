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

#ifndef ZF_LIB_CMN_TIME_H
#define ZF_LIB_CMN_TIME_H

/* -------------------------------------------------------------------------- */

/* return actual time */
long int time_actual(void);

/* set starting time and return its value */
long int time_start(void);
/* return time from starting point */
long int time_stop(void);

/* convert absolute time to formatted string */
void time_str(long, char*);
/* convert absolute time to formatted time string */
void time_str_time(long, char*);
/* convert absolute time to formatted date string */
void time_str_date(long, char*);
/* convert time interval in seconds to formatted string */
char* time_str_int(long);

/* -------------------------------------------------------------------------- */

#endif
