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

/* return internal ID of tethering potential
 
   s - the potential symbol */
short dlpoly_prm_teth_pot_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"harm")) return(DLPOLY_TETH_HM);
  if (str_compare(s,"rhrm")) return(DLPOLY_TETH_RH);
  if (str_compare(s,"quar")) return(DLPOLY_TETH_QR);
  return(id);
  }

/* convert internal ID of tethering potential to its symbol
 
   id - the internal ID
   s  - storage for the symbol */
char *dlpoly_prm_teth_pot_sym(short id) {
  static char s[80];
  switch (id) {
    case DLPOLY_TETH_HM: sprintf(s,"%s","harm"); break;
    case DLPOLY_TETH_RH: sprintf(s,"%s","rhrm"); break;
    case DLPOLY_TETH_QR: sprintf(s,"%s","quar"); break;
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */
