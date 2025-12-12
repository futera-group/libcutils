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
#include <cmn/message.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Return keyword name of topology bonding group

   id - internal ID of the bonding */
char *gmx_bond_key_name(short id) {
  static char key[80];
  switch (id) {
    case GMX_BTERM_BOND: sprintf(key,"%s","bonds");       break;
    case GMX_BTERM_ANGL: sprintf(key,"%s","angles");      break;
    case GMX_BTERM_DIHE: sprintf(key,"%s","dihedrals");   break;
    case GMX_BTERM_IMPR: sprintf(key,"%s","dihedrals");   break;
    case GMX_BTERM_PAIR: sprintf(key,"%s","pairs");       break;
    case GMX_BTERM_NBPR: sprintf(key,"%s","pairs");       break;
    case GMX_BTERM_CNST: sprintf(key,"%s","constraints"); break;
    case GMX_BTERM_STTL: sprintf(key,"%s","settles");     break;
    case GMX_BTERM_CMAP: sprintf(key,"%s","cmap");        break;
    case GMX_BTERM_EXCL: sprintf(key,"%s","exclusions");  break;
    default:
      msg_error_f("invalid bonding-group ID %d",1,id);
    }
  return(key);
  }

/* Return brief description of topology bonding group

   id - internal ID of the bonding */
char *gmx_bond_key_desc(short id) {
  static char key[80];
  switch (id) {
    case GMX_BTERM_BOND: sprintf(key,"%s","bond");            break;
    case GMX_BTERM_ANGL: sprintf(key,"%s","angle");           break;
    case GMX_BTERM_DIHE: sprintf(key,"%s","dihedral");        break;
    case GMX_BTERM_IMPR: sprintf(key,"%s","improper");        break;
    case GMX_BTERM_PAIR: sprintf(key,"%s","1-4 pair");        break;
    case GMX_BTERM_NBPR: sprintf(key,"%s","non-bonded pair"); break;
    case GMX_BTERM_CNST: sprintf(key,"%s","constraint");      break;
    case GMX_BTERM_STTL: sprintf(key,"%s","settle");          break;
    case GMX_BTERM_CMAP: sprintf(key,"%s","cmap correction"); break;
    case GMX_BTERM_EXCL: sprintf(key,"%s","exclusion");       break;
    default:
      msg_error_f("invalid bonding-group ID %d",1,id);
    }
  return(key);
  }

/* Return number of atoms involved in definition of specified bonding term

   id - internal ID of the bonding */
unsigned gmx_bond_key_natoms(short id) {
  switch (id) {
    case GMX_BTERM_BOND: return(2);
    case GMX_BTERM_ANGL: return(3);
    case GMX_BTERM_DIHE: return(4);
    case GMX_BTERM_IMPR: return(4);
    case GMX_BTERM_PAIR: return(2);
    case GMX_BTERM_NBPR: return(2);
    case GMX_BTERM_CNST: return(2);
    case GMX_BTERM_STTL: return(1);
    case GMX_BTERM_CMAP: return(5);
    }
  msg_error_f("invalid bonding-group ID %d",1,id);
  return(0);
  }

/* -------------------------------------------------------------------------- */
