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

/* return internal ID of four-body potential
 
   s - the potential symbol */
short dlpoly_prm_fbp_pot_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"harm")) return(DLPOLY_FBPT_HM);
  if (str_compare(s,"hcos")) return(DLPOLY_FBPT_HC);
  if (str_compare(s,"plan")) return(DLPOLY_FBPT_PL);
  return(id);
  }

/* convert internal ID of four-body potential to its symbol
 
   id - the internal ID
   s  - storage for the symbol */
char *dlpoly_prm_fbp_pot_sym(short id) {
  static char s[80];
  switch (id) {
    case DLPOLY_FBPT_HM: sprintf(s,"%s","harm"); break;
    case DLPOLY_FBPT_HC: sprintf(s,"%s","hcos"); break;
    case DLPOLY_FBPT_PL: sprintf(s,"%s","plan"); break;
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */
