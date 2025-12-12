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
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Return internal ID of particle type

   sym - particle-type symbol */
short gmx_ff_particle_type_id(char *sym) {
  if (str_compare(sym,"A"))
    return(GMX_PTYPE_ATOM);
  if (str_compare(sym,"S"))
    return(GMX_PTYPE_SHELL);
  if (str_compare(sym,"V") || str_compare(sym,"D"))
    return(GMX_PTYPE_VIRT);
  msg_error_f("unknown particle-type symbol \"%s\"",1,sym);
  return(0);
  }

/* Convert internal ID to particle-type symbol

   id - the internal particle-type ID */
char* gmx_ff_particle_type_name(short id) {
  static char name[80];
  name[0] = '\0';
  switch (id) {
    case GMX_PTYPE_ATOM:  sprintf(name,"%s","A"); break;
    case GMX_PTYPE_SHELL: sprintf(name,"%s","S"); break;
    case GMX_PTYPE_VIRT:  sprintf(name,"%s","V"); break;
    default:
      msg_error_f("unknown particle-type ID (%d)",1,id);
    }
  return(name);
  }

/* -------------------------------------------------------------------------- */
