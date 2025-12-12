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
#include <cmn/file.h>
#include <cmn/units.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* print out molecular coordinates of background charges to file

   d - pointer to gaussian data struct
   f - pointer to open output stream */
void gauss_mol_chrg_fprint(struct gauss_dat *d, FILE *f) {
  unsigned i;
  for (i=0; i<d->n_bg_charges; i++) {
    fprintf(f,"%12.6f %12.6f %12.6f %12.6f\n",
      d->bg_chrg[i][0]*CONV_B_ANG,
      d->bg_chrg[i][1]*CONV_B_ANG,
      d->bg_chrg[i][2]*CONV_B_ANG,
      d->bg_chrg[i][3]);
    }
  }

/* print out coordinates of background charges

   d - pointer to gaussian data struct
   f - name of the output file */
void gauss_mol_chrg_print(struct gauss_dat *d, char *f) {
  FILE *file = stdout;
  /* open file */
  if (f && f[0])
    file = file_open(f,"w");
  /* print charges */
  gauss_mol_chrg_fprint(d,file);
  /* close file */
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
