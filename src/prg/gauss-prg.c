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
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* return internal code of program version
 
   s - keyword from gaussian output file */
int gauss_prg_ver_id(char *s) {
  if (str_compare(s,"98"))
    return(GAUSS_PRG_98);
  else if (str_compare(s,"03"))
    return(GAUSS_PRG_03);
  else if (str_compare(s,"09"))
    return(GAUSS_PRG_09);
  return(GAUSS_PRG_UNK);
  }

/* return internal code of program revision
 
   s - keyword from gaussian output file */
int gauss_prg_rev_id(char *s) {
  if (str_compare(s,"A.01"))
    return(GAUSS_PRG_A01);
  else if (str_compare(s,"B.01"))
    return(GAUSS_PRG_B01);
  else if (str_compare(s,"C.01"))
    return(GAUSS_PRG_C01);
  return(GAUSS_PRG_UNK);
  }

/* -------------------------------------------------------------------------- */
