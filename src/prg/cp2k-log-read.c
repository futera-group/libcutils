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

#include <cmn/file.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* read data from cp2k output file
 
   d - pointer to cp2k data struct
   t - what to read
   f - name of the file */
void cp2k_log_read(struct cp2k_dat *d, int t, char *f) {
  FILE *file;
  file = file_open(f,"r");
  if (t & READ_HEADER)
    cp2k_log_read_header(d,file);
  if (t & READ_ATOMS)
    cp2k_log_read_atoms(d,file);
  if (t & READ_KINDS)
    cp2k_log_read_kinds(d,file);
  if (t & READ_OVERLAP)
    cp2k_log_read_overlap(d,file);
  if (t & READ_ENERGY)
    cp2k_log_read_energy(d,file);
  if (t & READ_CHARGES)
    cp2k_log_read_charges(d,file);
  if (t & READ_STATES)
    cp2k_log_read_states(d,file);
  if (t & READ_ECOUPL || t & READ_ECOUPL_INFO)
    cp2k_log_read_ecoupl(d,t,file);
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
