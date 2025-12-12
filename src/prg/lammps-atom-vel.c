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

#include <cmn/message.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* return number of velocity coefficients in specific atom style
 
   t - atom style ID */
unsigned lammps_atom_vel_nvar(short t) {
  unsigned n = 0;
  switch (t) {
    case LAMMPS_ATOM_ANGL: n = 4; break;
    case LAMMPS_ATOM_ATOM: n = 4; break;
    case LAMMPS_ATOM_BODY: n = 4; break;
    case LAMMPS_ATOM_BOND: n = 4; break;
    case LAMMPS_ATOM_CHRG: n = 4; break;
    case LAMMPS_ATOM_DIPL: n = 4; break;
    case LAMMPS_ATOM_DPDP: n = 4; break;
    case LAMMPS_ATOM_ELEC: n = 5; break;
    case LAMMPS_ATOM_ELLP: n = 7; break;
    case LAMMPS_ATOM_FULL: n = 4; break;
    case LAMMPS_ATOM_LINE: n = 4; break;
    case LAMMPS_ATOM_MESO: n = 4; break;
    case LAMMPS_ATOM_MOLE: n = 4; break;
    case LAMMPS_ATOM_PERI: n = 4; break;
    case LAMMPS_ATOM_SMDP: n = 4; break;
    case LAMMPS_ATOM_SPHR: n = 7; break;
    case LAMMPS_ATOM_TEMP: n = 4; break;
    case LAMMPS_ATOM_TRIP: n = 4; break;
    case LAMMPS_ATOM_WPCK: n = 4; break;
    default:
      msg_error_f("unknown atom style ID (%d)",1,t);
    }
  return(n);
  }

/* -------------------------------------------------------------------------- */
