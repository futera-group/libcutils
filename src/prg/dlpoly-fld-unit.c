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

/* return internal ID of force-field units
 
   s - the unit symbol */
short dlpoly_fld_unit_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"ev"))       return(DLPOLY_UNIT_EV);
  if (str_compare(s,"kcal"))     return(DLPOLY_UNIT_KC);
  if (str_compare(s,"kj"))       return(DLPOLY_UNIT_KJ);
  if (str_compare(s,"k"))        return(DLPOLY_UNIT_KB);
  if (str_compare(s,"internal")) return(DLPOLY_UNIT_IN);
  return(id);
  }

/* convert internal ID of force-field units to its symbol
 
   id - the internal ID */
char *dlpoly_fld_unit_sym(short id) {
  static char s[80];
  s[0] = '\0';
  /* full form keywords */
  switch (id) {
    case DLPOLY_UNIT_EV: sprintf(s,"%s","ev"); break;
    case DLPOLY_UNIT_KC: sprintf(s,"%s","kcal"); break;
    case DLPOLY_UNIT_KJ: sprintf(s,"%s","kj"); break;
    case DLPOLY_UNIT_KB: sprintf(s,"%s","k"); break;
    case DLPOLY_UNIT_IN: sprintf(s,"%s","internal"); break;
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */
