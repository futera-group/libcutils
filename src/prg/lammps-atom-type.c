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

/* return number of variables / coefficients in specific atom style
 
   t - atom style ID */
unsigned lammps_atom_type_nvar(short t) {
  unsigned n = 0;
  switch (t) {
    case LAMMPS_ATOM_ANGL: n =  6; break;
    case LAMMPS_ATOM_ATOM: n =  5; break;
    case LAMMPS_ATOM_BODY: n =  7; break;
    case LAMMPS_ATOM_BOND: n =  6; break;
    case LAMMPS_ATOM_CHRG: n =  6; break;
    case LAMMPS_ATOM_DIPL: n =  9; break;
    case LAMMPS_ATOM_DPDP: n =  6; break;
    case LAMMPS_ATOM_ELEC: n =  8; break;
    case LAMMPS_ATOM_ELLP: n =  7; break;
    case LAMMPS_ATOM_FULL: n =  7; break;
    case LAMMPS_ATOM_LINE: n =  8; break;
    case LAMMPS_ATOM_MESO: n =  8; break;
    case LAMMPS_ATOM_MOLE: n =  6; break;
    case LAMMPS_ATOM_PERI: n =  7; break;
    case LAMMPS_ATOM_SMDP: n = 10; break;
    case LAMMPS_ATOM_SPHR: n =  7; break;
    case LAMMPS_ATOM_TEMP: n =  8; break;
    case LAMMPS_ATOM_TRIP: n =  8; break;
    case LAMMPS_ATOM_WPCK: n = 11; break;
    default:
      msg_error_f("unknown atom style ID (%d)",1,t);
    }
  return(n);
  }

/* return internal ID of atom style
 
   s - the potential symbol */
short lammps_atom_type_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"angle"))      return(LAMMPS_ATOM_ANGL);
  if (str_compare(s,"atomic"))     return(LAMMPS_ATOM_ATOM);
  if (str_compare(s,"body"))       return(LAMMPS_ATOM_BODY);
  if (str_compare(s,"bond"))       return(LAMMPS_ATOM_BOND);
  if (str_compare(s,"charge"))     return(LAMMPS_ATOM_CHRG);
  if (str_compare(s,"dipole"))     return(LAMMPS_ATOM_DIPL);
  if (str_compare(s,"dpd"))        return(LAMMPS_ATOM_DPDP);
  if (str_compare(s,"electron"))   return(LAMMPS_ATOM_ELEC);
  if (str_compare(s,"ellipsoid"))  return(LAMMPS_ATOM_ELLP);
  if (str_compare(s,"full"))       return(LAMMPS_ATOM_FULL);
  if (str_compare(s,"line"))       return(LAMMPS_ATOM_LINE);
  if (str_compare(s,"meso"))       return(LAMMPS_ATOM_MESO);
  if (str_compare(s,"molecular"))  return(LAMMPS_ATOM_MOLE);
  if (str_compare(s,"peri"))       return(LAMMPS_ATOM_PERI);
  if (str_compare(s,"smd"))        return(LAMMPS_ATOM_SMDP);
  if (str_compare(s,"sphere"))     return(LAMMPS_ATOM_SPHR);
  if (str_compare(s,"template"))   return(LAMMPS_ATOM_TEMP);
  if (str_compare(s,"tri"))        return(LAMMPS_ATOM_TRIP);
  if (str_compare(s,"wavepacket")) return(LAMMPS_ATOM_WPCK);
  if (str_compare(s,"hybrid"))     return(LAMMPS_ATOM_HYBR);
  msg_error_f("unknown atom style \"%s\"",1,s);
  return(id);
  }

/* return symbol of atom specification type
 
   t - internal ID of the potential
   s - the potential symbol */
char* lammps_atom_type_sym(short t) {
  static char s[80];
  s[0]='\0';
  switch (t) {
    case LAMMPS_ATOM_ANGL: sprintf(s,"angle");      break;
    case LAMMPS_ATOM_ATOM: sprintf(s,"atomic");     break;
    case LAMMPS_ATOM_BODY: sprintf(s,"body");       break;
    case LAMMPS_ATOM_BOND: sprintf(s,"bond");       break;
    case LAMMPS_ATOM_CHRG: sprintf(s,"charge");     break;
    case LAMMPS_ATOM_DIPL: sprintf(s,"dipole");     break;
    case LAMMPS_ATOM_DPDP: sprintf(s,"dpd");        break;
    case LAMMPS_ATOM_ELEC: sprintf(s,"electron");   break;
    case LAMMPS_ATOM_ELLP: sprintf(s,"ellipsoid");  break;
    case LAMMPS_ATOM_FULL: sprintf(s,"full");       break;
    case LAMMPS_ATOM_LINE: sprintf(s,"line");       break;
    case LAMMPS_ATOM_MESO: sprintf(s,"meso");       break;
    case LAMMPS_ATOM_MOLE: sprintf(s,"molecular");  break;
    case LAMMPS_ATOM_PERI: sprintf(s,"peri");       break;
    case LAMMPS_ATOM_SMDP: sprintf(s,"smd");        break;
    case LAMMPS_ATOM_SPHR: sprintf(s,"sphere");     break;
    case LAMMPS_ATOM_TEMP: sprintf(s,"template");   break;
    case LAMMPS_ATOM_TRIP: sprintf(s,"tri");        break;
    case LAMMPS_ATOM_WPCK: sprintf(s,"wavepacket"); break;
    case LAMMPS_ATOM_HYBR: sprintf(s,"hybrid");     break;
    default:
      msg_error_f("unknown atom style (%d)",1,t);
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */
