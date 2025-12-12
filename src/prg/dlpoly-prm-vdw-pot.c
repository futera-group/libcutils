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

/* return internal ID of VdW potential
 
   s - the potential symbol */
short dlpoly_prm_vdw_pot_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"tab"))  return(DLPOLY_PAIR_TB);
  if (str_compare(s,"12-6")) return(DLPOLY_PAIR_12);
  if (str_compare(s,"lj"))   return(DLPOLY_PAIR_LJ);
  if (str_compare(s,"nm"))   return(DLPOLY_PAIR_NM);
  if (str_compare(s,"buck")) return(DLPOLY_PAIR_BC);
  if (str_compare(s,"bhm"))  return(DLPOLY_PAIR_BH);
  if (str_compare(s,"hbnd")) return(DLPOLY_PAIR_HB);
  if (str_compare(s,"snm"))  return(DLPOLY_PAIR_SF);
  if (str_compare(s,"mors")) return(DLPOLY_PAIR_MR);
  if (str_compare(s,"wca"))  return(DLPOLY_PAIR_WC);
  if (str_compare(s,"dpd"))  return(DLPOLY_PAIR_DP);
  if (str_compare(s,"amoe")) return(DLPOLY_PAIR_AM);
  return(id);
  }

/* convert internal ID of VdW potential to its symbol
 
   id - the internal ID
   s  - storage for the symbol */
char *dlpoly_prm_vdw_pot_sym(short id) {
  static char s[80];
  switch (id) {
    case DLPOLY_PAIR_TB: sprintf(s,"%s","tab");  break;
    case DLPOLY_PAIR_12: sprintf(s,"%s","12-6"); break;
    case DLPOLY_PAIR_LJ: sprintf(s,"%s","lj");   break;
    case DLPOLY_PAIR_NM: sprintf(s,"%s","nm");   break;
    case DLPOLY_PAIR_BC: sprintf(s,"%s","buck"); break;
    case DLPOLY_PAIR_BH: sprintf(s,"%s","bhm");  break;
    case DLPOLY_PAIR_HB: sprintf(s,"%s","hbnd"); break;
    case DLPOLY_PAIR_SF: sprintf(s,"%s","snm");  break;
    case DLPOLY_PAIR_MR: sprintf(s,"%s","mors"); break;
    case DLPOLY_PAIR_WC: sprintf(s,"%s","wca");  break;
    case DLPOLY_PAIR_DP: sprintf(s,"%s","dpd");  break;
    case DLPOLY_PAIR_AM: sprintf(s,"%s","amoe"); break;
    }
  return(s);
  }

/* return number of variables of VDW potential
 
   id - the potential specifier */
unsigned dlpoly_prm_vdw_pot_nvar(short id) {
  unsigned n = 0;
  switch (id) {
    case DLPOLY_PAIR_TB: n = 0; break;
    case DLPOLY_PAIR_12: n = 2; break;
    case DLPOLY_PAIR_LJ: n = 2; break;
    case DLPOLY_PAIR_NM: n = 4; break;
    case DLPOLY_PAIR_BC: n = 3; break;
    case DLPOLY_PAIR_BH: n = 5; break;
    case DLPOLY_PAIR_HB: n = 2; break;
    case DLPOLY_PAIR_SF: n = 5; break;
    case DLPOLY_PAIR_MR: n = 3; break;
    case DLPOLY_PAIR_WC: n = 3; break;
    case DLPOLY_PAIR_DP: n = 2; break;
    case DLPOLY_PAIR_AM: n = 2; break;
    }
  return(n);
  }

/* -------------------------------------------------------------------------- */
