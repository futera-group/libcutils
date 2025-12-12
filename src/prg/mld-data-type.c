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
#include "prg/molden.h"

/* -------------------------------------------------------------------------- */

/* convert data type symbol to internal ID

   s - symbol of the data type */
short mld_data_type_id(char *s) {
  if (str_compare(s,"Atoms"))         return(MLD_DATA_ATOM);
  if (str_compare(s,"GTO"))           return(MLD_DATA_GTO);
  if (str_compare(s,"STO"))           return(MLD_DATA_STO);
  if (str_compare(s,"MO"))            return(MLD_DATA_MO);
  if (str_compare(s,"5D"))            return(MLD_DATA_5D);
  if (str_compare(s,"5D7F"))          return(MLD_DATA_5D7F);
  if (str_compare(s,"5D10F"))         return(MLD_DATA_5D10F);
  if (str_compare(s,"7F"))            return(MLD_DATA_7F);
  if (str_compare(s,"9G"))            return(MLD_DATA_9G);
  if (str_compare(s,"SCFCONV"))       return(MLD_DATA_SCF);
  if (str_compare(s,"GEOCONV"))       return(MLD_DATA_GEO);
  if (str_compare(s,"GEOMETRIES"))    return(MLD_DATA_GEOC);
  if (str_compare(s,"FREQ"))          return(MLD_DATA_FREQ);
  if (str_compare(s,"FR-COORD"))      return(MLD_DATA_FRC);
  if (str_compare(s,"FR-NORM-COORD")) return(MLD_DATA_FRNC);
  if (str_compare(s,"INT"))           return(MLD_DATA_FRI);
  return(0);
  }

/* convert data type internal ID to symbol

   t - internal code of the data type */
char* mld_data_type_name(short t) {
  static char name[80];
  switch (t) {
    case MLD_DATA_ATOM:  sprintf(name,"Atoms");         break;
    case MLD_DATA_GTO:   sprintf(name,"GTO");           break;
    case MLD_DATA_STO:   sprintf(name,"STO");           break;
    case MLD_DATA_MO:    sprintf(name,"MO");            break;
    case MLD_DATA_5D:    sprintf(name,"5D");            break;
    case MLD_DATA_5D7F:  sprintf(name,"5D7F");          break;
    case MLD_DATA_5D10F: sprintf(name,"5D10F");         break;
    case MLD_DATA_7F:    sprintf(name,"7F");            break;
    case MLD_DATA_9G:    sprintf(name,"9G");            break;
    case MLD_DATA_SCF:   sprintf(name,"SCFCONV");       break;
    case MLD_DATA_GEO:   sprintf(name,"GEOCONV");       break;
    case MLD_DATA_GEOC:  sprintf(name,"GEOMETRIES");    break;
    case MLD_DATA_FREQ:  sprintf(name,"FREQ");          break;
    case MLD_DATA_FRC:   sprintf(name,"FR-COORD");      break;
    case MLD_DATA_FRNC:  sprintf(name,"FR-NORM-COORD"); break;
    case MLD_DATA_FRI:   sprintf(name,"INT");           break;
    default: sprintf(name,"X");                         break;
    }
  return(name);
  }

/* -------------------------------------------------------------------------- */
