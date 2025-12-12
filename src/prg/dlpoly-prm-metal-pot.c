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

/* return internal ID of metal potential
 
   s - the potential symbol */
short dlpoly_prm_metal_pot_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"eam"))  return(DLPOLY_METL_EA);
  if (str_compare(s,"eeam")) return(DLPOLY_METL_EE);
  if (str_compare(s,"2bea")) return(DLPOLY_METL_2B);
  if (str_compare(s,"2bee")) return(DLPOLY_METL_2E);
  if (str_compare(s,"fnsc")) return(DLPOLY_METL_FS);
  if (str_compare(s,"exfs")) return(DLPOLY_METL_EF);
  if (str_compare(s,"stch")) return(DLPOLY_METL_SC);
  if (str_compare(s,"gupt")) return(DLPOLY_METL_GP);
  if (str_compare(s,"mbpc")) return(DLPOLY_METL_MB);
  return(id);
  }

/* convert internal ID of metal potential to its symbol
 
   id - the internal ID
   s  - storage for the symbol */
char *dlpoly_prm_metal_pot_sym(short id) {
  static char s[80];
  switch (id) {
    case DLPOLY_METL_EA: sprintf(s,"%s","eam");  break;
    case DLPOLY_METL_EE: sprintf(s,"%s","eeam"); break;
    case DLPOLY_METL_2B: sprintf(s,"%s","2bea"); break;
    case DLPOLY_METL_2E: sprintf(s,"%s","2bee"); break;
    case DLPOLY_METL_FS: sprintf(s,"%s","fnsc"); break;
    case DLPOLY_METL_EF: sprintf(s,"%s","exfs"); break;
    case DLPOLY_METL_SC: sprintf(s,"%s","stch"); break;
    case DLPOLY_METL_GP: sprintf(s,"%s","gupt"); break;
    case DLPOLY_METL_MB: sprintf(s,"%s","mbpc"); break;
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */
