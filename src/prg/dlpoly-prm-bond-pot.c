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

/* return internal ID of bond potential
 
   s - the potential symbol */
short dlpoly_prm_bond_pot_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"harm")) return(DLPOLY_BOND_HM);
  if (str_compare(s,"-hrm")) return(DLPOLY_BOND_HM);
  if (str_compare(s,"mors")) return(DLPOLY_BOND_MR);
  if (str_compare(s,"-mrs")) return(DLPOLY_BOND_MR);
  if (str_compare(s,"12-6")) return(DLPOLY_BOND_12);
  if (str_compare(s,"-126")) return(DLPOLY_BOND_12);
  if (str_compare(s,"lj"))   return(DLPOLY_BOND_LJ);
  if (str_compare(s,"-lj"))  return(DLPOLY_BOND_LJ);
  if (str_compare(s,"quar")) return(DLPOLY_BOND_QR);
  if (str_compare(s,"-qur")) return(DLPOLY_BOND_QR);
  if (str_compare(s,"buck")) return(DLPOLY_BOND_BC);
  if (str_compare(s,"-bck")) return(DLPOLY_BOND_BC);
  if (str_compare(s,"coul")) return(DLPOLY_BOND_CL);
  if (str_compare(s,"-cul")) return(DLPOLY_BOND_CL);
  if (str_compare(s,"fene")) return(DLPOLY_BOND_SF);
  if (str_compare(s,"-fne")) return(DLPOLY_BOND_SF);
  if (str_compare(s,"amoe")) return(DLPOLY_BOND_AM);
  if (str_compare(s,"-amo")) return(DLPOLY_BOND_AM);
  return(id);
  }

/* convert internal ID of bond potential to its symbol
 
   id - the internal ID
   m  - minus form of the keyword
   s  - storage for the symbol */
char *dlpoly_prm_bond_pot_sym(short id, short m) {
  static char s[80];
  if (m) {
    /* full form keywords */
    switch (id) {
      case DLPOLY_BOND_HM: sprintf(s,"%s","-hrm"); break;
      case DLPOLY_BOND_MR: sprintf(s,"%s","-mrs"); break;
      case DLPOLY_BOND_12: sprintf(s,"%s","-126"); break;
      case DLPOLY_BOND_LJ: sprintf(s,"%s","-lj");  break;
      case DLPOLY_BOND_RH: sprintf(s,"%s","-rhm"); break;
      case DLPOLY_BOND_QR: sprintf(s,"%s","-qur"); break;
      case DLPOLY_BOND_BC: sprintf(s,"%s","-bck"); break;
      case DLPOLY_BOND_CL: sprintf(s,"%s","-cul"); break;
      case DLPOLY_BOND_SF: sprintf(s,"%s","-fne"); break;
      case DLPOLY_BOND_AM: sprintf(s,"%s","-amo"); break;
      }
    }
  else {
    /* minus form keywords */
    switch (id) {
      case DLPOLY_BOND_HM: sprintf(s,"%s","harm"); break;
      case DLPOLY_BOND_MR: sprintf(s,"%s","mors"); break;
      case DLPOLY_BOND_12: sprintf(s,"%s","12-6"); break;
      case DLPOLY_BOND_LJ: sprintf(s,"%s","lj");   break;
      case DLPOLY_BOND_RH: sprintf(s,"%s","rhrm"); break;
      case DLPOLY_BOND_QR: sprintf(s,"%s","quar"); break;
      case DLPOLY_BOND_BC: sprintf(s,"%s","buck"); break;
      case DLPOLY_BOND_CL: sprintf(s,"%s","coul"); break;
      case DLPOLY_BOND_SF: sprintf(s,"%s","fene"); break;
      case DLPOLY_BOND_AM: sprintf(s,"%s","amoe"); break;
      }
    }
  return(s);
  }

/* return number of variables of specific bond potential
 
   id - the potential specifier */
unsigned dlpoly_prm_bond_pot_nvar(short id) {
  unsigned n = 0;
  switch (id) {
    case DLPOLY_BOND_HM: n = 2; break;
    case DLPOLY_BOND_MR: n = 3; break;
    case DLPOLY_BOND_12: n = 2; break;
    case DLPOLY_BOND_LJ: n = 2; break;
    case DLPOLY_BOND_RH: n = 3; break;
    case DLPOLY_BOND_QR: n = 4; break;
    case DLPOLY_BOND_BC: n = 3; break;
    case DLPOLY_BOND_CL: n = 1; break;
    case DLPOLY_BOND_SF: n = 3; break;
    case DLPOLY_BOND_AM: n = 2; break;
    }
  return(n);
  }

/* -------------------------------------------------------------------------- */
