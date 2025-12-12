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
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* return internal ID of tersoff potential
 
   s - the potential symbol */
short dlpoly_prm_ters_pot_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"ters")) return(DLPOLY_TERS_TR);
  if (str_compare(s,"kihs")) return(DLPOLY_TERS_KH);
  return(id);
  }

/* convert internal ID of tersoff potential to its symbol
 
   id - the internal ID
   s  - storage for the symbol */
char *dlpoly_prm_ters_pot_sym(short id) {
  static char s[80];
  switch (id) {
    case DLPOLY_TERS_TR: sprintf(s,"%s","ters"); break;
    case DLPOLY_TERS_KH: sprintf(s,"%s","kihs"); break;
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */
