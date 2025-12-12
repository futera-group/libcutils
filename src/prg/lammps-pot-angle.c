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

/* return number of variables / coefficients in specific angle potential
 
   t - type of the potential */
unsigned lammps_pot_angle_nvar(short t) {
  unsigned n = 0;
  switch (t) {
    case LAMMPS_ANGL_CHRM: n = 4; break;
    case LAMMPS_ANGL_CLS2: n = 4; break;
    case LAMMPS_ANGL_COSN: n = 1; break;
    case LAMMPS_ANGL_COSD: n = 2; break;
    case LAMMPS_ANGL_COSP: n = 3; break;
    case LAMMPS_ANGL_COSS: n = 2; break;
    case LAMMPS_ANGL_HARM: n = 2; break;
    default:
      msg_error_f("unknown angle potential type (%d)",1,t);
    }
  return(n);
  }

/* return internal ID of angle potential
 
   s - the potential symbol */
short lammps_pot_angle_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"none"))            return(LAMMPS_ANGL_NONE);
  if (str_compare(s,"zero"))            return(LAMMPS_ANGL_ZERO);
  if (str_compare(s,"hybrid"))          return(LAMMPS_ANGL_HYBR);
  if (str_compare(s,"charmm"))          return(LAMMPS_ANGL_CHRM);
  if (str_compare(s,"class2"))          return(LAMMPS_ANGL_CLS2);
  if (str_compare(s,"cosine"))          return(LAMMPS_ANGL_COSN);
  if (str_compare(s,"cosine/delta"))    return(LAMMPS_ANGL_COSD);
  if (str_compare(s,"cosine/periodic")) return(LAMMPS_ANGL_COSP);
  if (str_compare(s,"cosine/squared"))  return(LAMMPS_ANGL_COSS);
  if (str_compare(s,"harmonic"))        return(LAMMPS_ANGL_HARM);
  if (str_compare(s,"table"))           return(LAMMPS_ANGL_TABL);
  msg_error_f("unknown angle potential type \"%s\"",1,s);
  return(id);
  }

/* return symbol of angle potential
 
   t - internal ID of the potential
   s - the potential symbol */
char* lammps_pot_angle_sym(short t) {
  static char s[80];
  s[0]='\0';
  switch (t) {
    case LAMMPS_ANGL_NONE: sprintf(s,"none");            break;
    case LAMMPS_ANGL_ZERO: sprintf(s,"zero");            break;
    case LAMMPS_ANGL_HYBR: sprintf(s,"hybrid");          break;
    case LAMMPS_ANGL_CHRM: sprintf(s,"charmm");          break;
    case LAMMPS_ANGL_CLS2: sprintf(s,"class2");          break;
    case LAMMPS_ANGL_COSN: sprintf(s,"cosine");          break;
    case LAMMPS_ANGL_COSD: sprintf(s,"cosine/delta");    break;
    case LAMMPS_ANGL_COSP: sprintf(s,"cosine/periodic"); break;
    case LAMMPS_ANGL_COSS: sprintf(s,"cosine/squared");  break;
    case LAMMPS_ANGL_HARM: sprintf(s,"harmonic");        break;
    case LAMMPS_ANGL_TABL: sprintf(s,"table");           break;
    default:
      msg_error_f("unknown angle potential type (%d)",1,t);
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */
