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

#include "qmc/gto.h"

/* -------------------------------------------------------------------------- */

/* get name of the 2e integral class

   a,b,c,d - angular momenta of (ab|cd) integral */
char *gto_int_2e_name(unsigned a, unsigned b, unsigned c, unsigned d) {
  static char name[10] = "(XX|XX)";
  name[1] = gto_ang_name(a);
  name[2] = gto_ang_name(b);
  name[4] = gto_ang_name(c);
  name[5] = gto_ang_name(d);
  return(name);
  }

/* -------------------------------------------------------------------------- */
