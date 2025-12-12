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

#include <cmn/string.h>
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* return internal ID of dihedral angle potential
 
   s - the potential symbol */
short dlpoly_prm_dihed_pot_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"cos"))  return(DLPOLY_DIHE_CS);
  if (str_compare(s,"harm")) return(DLPOLY_DIHE_HM);
  if (str_compare(s,"hcos")) return(DLPOLY_DIHE_HC);
  if (str_compare(s,"cos3")) return(DLPOLY_DIHE_TC);
  if (str_compare(s,"ryck")) return(DLPOLY_DIHE_RB);
  if (str_compare(s,"rbf"))  return(DLPOLY_DIHE_FR);
  if (str_compare(s,"opls")) return(DLPOLY_DIHE_OP);
  return(id);
  }

/* convert internal ID of dihedral angle potential to its symbol
 
   id - the internal ID
   s  - storage for the symbol */
char *dlpoly_prm_dihed_pot_sym(short id) {
  static char s[80];
  switch (id) {
    case DLPOLY_DIHE_CS: sprintf(s,"%s","cos");  break;
    case DLPOLY_DIHE_HM: sprintf(s,"%s","harm"); break;
    case DLPOLY_DIHE_HC: sprintf(s,"%s","hcos"); break;
    case DLPOLY_DIHE_TC: sprintf(s,"%s","cos3"); break;
    case DLPOLY_DIHE_RB: sprintf(s,"%s","ryck"); break;
    case DLPOLY_DIHE_FR: sprintf(s,"%s","rbf");  break;
    case DLPOLY_DIHE_OP: sprintf(s,"%s","opls"); break;
    }
  return(s);
  }

/* return number of variables of specific dihedral angle potential
 
   id - the potential specifier */
unsigned dlpoly_prm_dihed_pot_nvar(short id) {
  unsigned n = 0;
  switch (id) {
    case DLPOLY_DIHE_CS: n = 3; break;
    case DLPOLY_DIHE_HM: n = 2; break;
    case DLPOLY_DIHE_HC: n = 2; break;
    case DLPOLY_DIHE_TC: n = 3; break;
    case DLPOLY_DIHE_RB: n = 1; break;
    case DLPOLY_DIHE_FR: n = 1; break;
    case DLPOLY_DIHE_OP: n = 3; break;
    }
  return(n);
  }

/* -------------------------------------------------------------------------- */
