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
#include <mol/atom.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* print out molecular structure coordinates to file

   d - pointer to gaussian data struct
   f - pointer to open output stream */
void gauss_mol_coord_fprint(struct gauss_dat *d, FILE *f) {
  unsigned i;
  for (i=0; i<d->n_atoms; i++) {
    fprintf(f,"%5s %12.6f %12.6f %12.6f\n",atom_name(d->atom[i].num),
      d->atom[i].coord[0]*CONV_B_ANG,
      d->atom[i].coord[1]*CONV_B_ANG,
      d->atom[i].coord[2]*CONV_B_ANG);
    }
  }

/* print out molecular structure coordinates

   d - pointer to gaussian data struct
   f - name of the output file */
void gauss_mol_coord_print(struct gauss_dat *d, char *f) {
  FILE *file = stdout;
  /* open file */
  if (f && f[0])
    file = file_open(f,"w");
  /* print coordinates */
  gauss_mol_coord_fprint(d,file);
  /* close file */
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
