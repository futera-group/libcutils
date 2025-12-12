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

/* return number of variables / coefficients in specific dihedral potential
 
   t - type of the potential */
unsigned lammps_pot_dihed_nvar(short t) {
  unsigned n = 0;
  switch (t) {
    case LAMMPS_DIHE_CHRM: n = 4; break;
    case LAMMPS_DIHE_CLS2: n = 6; break;
    case LAMMPS_DIHE_HARM: n = 3; break;
    case LAMMPS_DIHE_HELX: n = 3; break;
    case LAMMPS_DIHE_MULT: n = 5; break;
    case LAMMPS_DIHE_OPLS: n = 4; break;
    default:
      msg_error_f("unknown dihedral-angle potential type (%d)",1,t);
    }
  return(n);
  }

/* return internal ID of dihedral-angle potential
 
   s - the potential symbol */
short lammps_pot_dihed_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"none"))           return(LAMMPS_DIHE_NONE);
  if (str_compare(s,"zero"))           return(LAMMPS_DIHE_ZERO);
  if (str_compare(s,"hybrid"))         return(LAMMPS_DIHE_HYBR);
  if (str_compare(s,"charmm"))         return(LAMMPS_DIHE_CHRM);
  if (str_compare(s,"class2"))         return(LAMMPS_DIHE_CLS2);
  if (str_compare(s,"harmonic"))       return(LAMMPS_DIHE_HARM);
  if (str_compare(s,"helix"))          return(LAMMPS_DIHE_HELX);
  if (str_compare(s,"multi/harmonic")) return(LAMMPS_DIHE_MULT);
  if (str_compare(s,"opls"))           return(LAMMPS_DIHE_OPLS);
  msg_error_f("unknown dihedral-angle potential type \"%s\"",1,s);
  return(id);
  }

/* return symbol of dihedral-angle potential
 
   t - internal ID of the potential
   s - the potential symbol */
char* lammps_pot_dihed_sym(short t) {
  static char s[80];
  s[0]='\0';
  switch (t) {
    case LAMMPS_DIHE_NONE: sprintf(s,"none");           break;
    case LAMMPS_DIHE_ZERO: sprintf(s,"zero");           break;
    case LAMMPS_DIHE_HYBR: sprintf(s,"hybrid");         break;
    case LAMMPS_DIHE_CHRM: sprintf(s,"charmm");         break;
    case LAMMPS_DIHE_CLS2: sprintf(s,"class2");         break;
    case LAMMPS_DIHE_HARM: sprintf(s,"harmonic");       break;
    case LAMMPS_DIHE_HELX: sprintf(s,"helix");          break;
    case LAMMPS_DIHE_MULT: sprintf(s,"multi/harmonic"); break;
    case LAMMPS_DIHE_OPLS: sprintf(s,"opls");           break;
    default:
      msg_error_f("unknown dihedral-angle potential type (%d)",1,t);
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */
