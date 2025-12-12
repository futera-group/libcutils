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

/* return internal ID of three-body potential
 
   s - the potential symbol */
short dlpoly_prm_tbp_pot_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"harm")) return(DLPOLY_TBPT_HM);
  if (str_compare(s,"thrm")) return(DLPOLY_TBPT_TH);
  if (str_compare(s,"shrm")) return(DLPOLY_TBPT_SH);
  if (str_compare(s,"bvs1")) return(DLPOLY_TBPT_SV);
  if (str_compare(s,"bvs2")) return(DLPOLY_TBPT_TV);
  if (str_compare(s,"hbnd")) return(DLPOLY_TBPT_HB);
  return(id);
  }

/* convert internal ID of three-body potential to its symbol
 
   id - the internal ID
   s  - storage for the symbol */
char *dlpoly_prm_tbp_pot_sym(short id) {
  static char s[80];
  switch (id) {
    case DLPOLY_TBPT_HM: sprintf(s,"%s","harm"); break;
    case DLPOLY_TBPT_TH: sprintf(s,"%s","thrm"); break;
    case DLPOLY_TBPT_SH: sprintf(s,"%s","shrm"); break;
    case DLPOLY_TBPT_SV: sprintf(s,"%s","bvs1"); break;
    case DLPOLY_TBPT_TV: sprintf(s,"%s","bvs2"); break;
    case DLPOLY_TBPT_HB: sprintf(s,"%s","hbnd"); break;
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */
