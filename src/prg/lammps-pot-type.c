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
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* return internal ID of potential type
 
   s - the potential symbol */
short lammps_pot_type_id(char *s) {
  short id = 0;
  str_lowcase(s);
  if (str_compare(s,"pair"))              return(LAMMPS_POT_PAIR);
  if (str_compare(s,"pairij"))            return(LAMMPS_POT_PRIJ);
  if (str_compare(s,"bond"))              return(LAMMPS_POT_BOND);
  if (str_compare(s,"angle"))             return(LAMMPS_POT_ANGL);
  if (str_compare(s,"dihedral"))          return(LAMMPS_POT_DIHE);
  if (str_compare(s,"improper"))          return(LAMMPS_POT_IMPR);
  if (str_compare(s,"bondbond"))          return(LAMMPS_POT_BNBN);
  if (str_compare(s,"bondangle"))         return(LAMMPS_POT_BNAN);
  if (str_compare(s,"middlebondtorsion")) return(LAMMPS_POT_MBTR);
  if (str_compare(s,"endbondtorsion"))    return(LAMMPS_POT_EBTR);
  if (str_compare(s,"angletorsion"))      return(LAMMPS_POT_ANTR);
  if (str_compare(s,"angleangle"))        return(LAMMPS_POT_ANAN);
  if (str_compare(s,"angleangletorsion")) return(LAMMPS_POT_AATR);
  if (str_compare(s,"bondbond13"))        return(LAMMPS_POT_BB13);
  return(id);
  }

/* -------------------------------------------------------------------------- */
