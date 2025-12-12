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

/* return internal ID of external-field potential
 
   s - the potential symbol */
short dlpoly_prm_efld_pot_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"elec")) return(DLPOLY_EFLD_EL);
  if (str_compare(s,"oshr")) return(DLPOLY_EFLD_OS);
  if (str_compare(s,"shrx")) return(DLPOLY_EFLD_CS);
  if (str_compare(s,"grav")) return(DLPOLY_EFLD_GV);
  if (str_compare(s,"magn")) return(DLPOLY_EFLD_MG);
  if (str_compare(s,"sphr")) return(DLPOLY_EFLD_SP);
  if (str_compare(s,"zbnd")) return(DLPOLY_EFLD_RW);
  if (str_compare(s,"xpis")) return(DLPOLY_EFLD_XP);
  if (str_compare(s,"zres")) return(DLPOLY_EFLD_HA);
  if (str_compare(s,"zrs-")) return(DLPOLY_EFLD_HB);
  if (str_compare(s,"zrs+")) return(DLPOLY_EFLD_HC);
  if (str_compare(s,"osel")) return(DLPOLY_EFLD_OE);
  return(id);
  }

/* convert internal ID of external-field potential to its symbol
 
   id - the internal ID
   s  - storage for the symbol */
char *dlpoly_prm_efld_pot_sym(short id) {
  static char s[80];
  switch (id) {
    case DLPOLY_EFLD_EL: sprintf(s,"%s","elec"); break;
    case DLPOLY_EFLD_OS: sprintf(s,"%s","oshr"); break;
    case DLPOLY_EFLD_CS: sprintf(s,"%s","shrx"); break;
    case DLPOLY_EFLD_GV: sprintf(s,"%s","grav"); break;
    case DLPOLY_EFLD_MG: sprintf(s,"%s","magn"); break;
    case DLPOLY_EFLD_SP: sprintf(s,"%s","sphr"); break;
    case DLPOLY_EFLD_RW: sprintf(s,"%s","zbnd"); break;
    case DLPOLY_EFLD_XP: sprintf(s,"%s","xpis"); break;
    case DLPOLY_EFLD_HA: sprintf(s,"%s","zres"); break;
    case DLPOLY_EFLD_HB: sprintf(s,"%s","zrs-"); break;
    case DLPOLY_EFLD_HC: sprintf(s,"%s","zrs+"); break;
    case DLPOLY_EFLD_OE: sprintf(s,"%s","osel"); break;
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */
