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
#include <cmn/string.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* return number of variables / coefficients in specific bond potential
 
   t - type of the potential */
unsigned lammps_pot_bond_nvar(short t) {
  unsigned n = 0;
  switch (t) {
    case LAMMPS_BOND_CLS2: n = 4; break;
    case LAMMPS_BOND_FENE: n = 4; break;
    case LAMMPS_BOND_FEEX: n = 4; break;
    case LAMMPS_BOND_HARM: n = 2; break;
    case LAMMPS_BOND_MORS: n = 3; break;
    case LAMMPS_BOND_NONL: n = 3; break;
    case LAMMPS_BOND_QUAR: n = 5; break;
    default:
      msg_error_f("unknown bond potential type (%d)",1,t);
    }
  return(n);
  }

/* return internal ID of bond potential
 
   s - the potential symbol */
short lammps_pot_bond_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"none"))        return(LAMMPS_BOND_NONE);
  if (str_compare(s,"hybrid"))      return(LAMMPS_BOND_HYBR);
  if (str_compare(s,"class2"))      return(LAMMPS_BOND_CLS2);
  if (str_compare(s,"fene"))        return(LAMMPS_BOND_FENE);
  if (str_compare(s,"fene/expand")) return(LAMMPS_BOND_FEEX);
  if (str_compare(s,"harmonic"))    return(LAMMPS_BOND_HARM);
  if (str_compare(s,"morse"))       return(LAMMPS_BOND_MORS);
  if (str_compare(s,"nonlinear"))   return(LAMMPS_BOND_NONL);
  if (str_compare(s,"quartic"))     return(LAMMPS_BOND_QUAR);
  if (str_compare(s,"table"))       return(LAMMPS_BOND_TABL);
  msg_error_f("unknown bond potential type \"%s\"",1,s);
  return(id);
  }

/* return symbol of bond potential
 
   t - internal ID of the potential
   s - the potential symbol */
char* lammps_pot_bond_sym(short t) {
  static char s[80];
  s[0]='\0';
  switch (t) {
    case LAMMPS_BOND_NONE: sprintf(s,"none");        break;
    case LAMMPS_BOND_HYBR: sprintf(s,"hybrid");      break;
    case LAMMPS_BOND_CLS2: sprintf(s,"class2");      break;
    case LAMMPS_BOND_FENE: sprintf(s,"fene");        break;
    case LAMMPS_BOND_FEEX: sprintf(s,"fene/expand"); break;
    case LAMMPS_BOND_HARM: sprintf(s,"harmonic");    break;
    case LAMMPS_BOND_MORS: sprintf(s,"morse");       break;
    case LAMMPS_BOND_NONL: sprintf(s,"nonlinear");   break;
    case LAMMPS_BOND_QUAR: sprintf(s,"quartic");     break;
    case LAMMPS_BOND_TABL: sprintf(s,"table");       break;
    default:
      msg_error_f("unknown bond potential type (%d)",1,t);
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */
