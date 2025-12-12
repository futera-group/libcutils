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
#include <cmn/string.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* convert calculation type symbol to internal ID

   s - symbol of calculation type */
short gauss_job_type_id(char *s) {
  if (str_compare(s,"SP"))                   return(GAUSS_JOB_SP);
  if (str_compare(s,"FOPT"))                 return(GAUSS_JOB_FOPT);
  if (str_compare(s,"POPT"))                 return(GAUSS_JOB_POPT);
  if (str_compare(s,"FTS"))                  return(GAUSS_JOB_FTS);
  if (str_compare(s,"PTS"))                  return(GAUSS_JOB_PTS);
  if (str_compare(s,"FSADDLE"))              return(GAUSS_JOB_FSDD);
  if (str_compare(s,"PSADDLE"))              return(GAUSS_JOB_PSDD);
  if (str_compare(s,"FORCE"))                return(GAUSS_JOB_FORCE);
  if (str_compare(s,"FREQ"))                 return(GAUSS_JOB_FREQ);
  if (str_compare(s,"SCAN"))                 return(GAUSS_JOB_SCAN);
  if (str_compare(s,"GUESS=ONLY"))           return(GAUSS_JOB_GUESS);
  if (str_compare(s,"LST"))                  return(GAUSS_JOB_LST);
  if (str_compare(s,"STABILITY"))            return(GAUSS_JOB_STAB);
  if (str_compare(s,"REARCHIVE/MS-RESTART")) return(GAUSS_JOB_RARCH);
  if (str_compare(s,"MIXED"))                return(GAUSS_JOB_MIXED);
  return(0);
  }

/* convert calculation type internal ID to symbol

   s - symbol of calculation type (output)
   t - internal code of calculation type */
char* gauss_job_type_name(short t) {
  static char name[80];
  switch (t) {
    case GAUSS_JOB_SP:    sprintf(name,"SP"); break;
    case GAUSS_JOB_FOPT:  sprintf(name,"FOPT"); break;
    case GAUSS_JOB_POPT:  sprintf(name,"POPT"); break;
    case GAUSS_JOB_FTS:   sprintf(name,"FTS"); break;
    case GAUSS_JOB_PTS:   sprintf(name,"PTS"); break;
    case GAUSS_JOB_FSDD:  sprintf(name,"FSADDLE"); break;
    case GAUSS_JOB_PSDD:  sprintf(name,"PSADDLE"); break;
    case GAUSS_JOB_FORCE: sprintf(name,"FORCE"); break;
    case GAUSS_JOB_FREQ:  sprintf(name,"FREQ"); break;
    case GAUSS_JOB_SCAN:  sprintf(name,"SCAN"); break;
    case GAUSS_JOB_GUESS: sprintf(name,"GUESS=ONLY"); break;
    case GAUSS_JOB_LST:   sprintf(name,"LST"); break;
    case GAUSS_JOB_STAB:  sprintf(name,"STABILITY"); break;
    case GAUSS_JOB_RARCH: sprintf(name,"REARCHIVE/MS-RESTART"); break;
    case GAUSS_JOB_MIXED: sprintf(name,"MIXED"); break;
    default: sprintf(name,"X"); break;
    }
  return(name);
  }

/* -------------------------------------------------------------------------- */
