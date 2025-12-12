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
#include <cmn/message.h>
#include <cmn/string.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* return internal ID of output type 
 
   key - print level keyword */
short cp2k_print_level_id(char *key) {
  char s[80] = "\n";
  str_lowcase_copy(key,s);
  if (str_compare(s,"debug"))
    return(OUT_TYPE_DEBG);
  if (str_compare(s,"high"))
    return(OUT_TYPE_HIGH);
  if (str_compare(s,"low"))
    return(OUT_TYPE_LOWX);
  if (str_compare(s,"medium"))
    return(OUT_TYPE_MEDM);
  if (str_compare(s,"silent"))
    return(OUT_TYPE_SLNT);
  msg_error_f("unknown output type keyword \"%s\"",1,key);
  return(0);
  }

/* convert internal output type ID to corresponding keyword
 
   key - storage for the keyword
   id  - the internal ID */
void cp2k_print_level_name(char *key, short id) {
  switch (id) {
    case OUT_TYPE_DEBG: sprintf(key,"debug");  break;
    case OUT_TYPE_HIGH: sprintf(key,"high");   break;
    case OUT_TYPE_LOWX: sprintf(key,"low");    break;
    case OUT_TYPE_MEDM: sprintf(key,"medium"); break;
    case OUT_TYPE_SLNT: sprintf(key,"silent"); break;
    default:
      msg_error_f("unknown internal output-type ID (%d)",1,id);
    }                      
  }

/* -------------------------------------------------------------------------- */
