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

/* write vibrational-frequency normal modes to molden-type file
 
   d    - gaussian data
   file - name of the file */
void gauss_freq_write(struct gauss_dat *d, char *file) {
  unsigned i,j,k;
  FILE *f;
  /* open data file */
  f = file_open(file,"w");
  /* header */
  fprintf(f,"[Molden Format]\n");
  /* coordinates [au] */
  fprintf(f,"[Atoms] Angs\n");
  for (i=0; i<d->n_atoms; i++) {
    fprintf(f," %-3s %5d %3d",atom_name(d->atom[i].num),i+1,d->atom[i].num);
    for (j=0; j<3; j++)
      fprintf(f," %20.10E",d->atom[i].coord[j]*CONV_B_ANG);
    fprintf(f,"\n");
    }
  /* frequencies [cm-1] */
  fprintf(f,"[FREQ]\n");
  for (i=0; i<d->n_vib_modes; i++)
    fprintf(f,"%20.10E\n",d->freq[i].freq);
  /* coordinates [au] */
  fprintf(f,"[FR-COORD]\n");
  for (i=0; i<d->n_atoms; i++) {
    fprintf(f," %-3s",atom_name(d->atom[i].num));
    for (j=0; j<3; j++)
      fprintf(f," %20.10E",d->atom[i].coord[j]);
    fprintf(f,"\n");
    }
  /* normal modes [au] */
  fprintf(f,"[FR-NORM-COORD]\n");
  for (i=0; i<d->n_vib_modes; i++) {
    fprintf(f,"vibration %d\n",i+1);
    for (j=0; j<d->n_atoms; j++) {
      for (k=0; k<3; k++)
        fprintf(f," %20.10E",d->freq[i].displ[j][k]);
      fprintf(f,"\n");
      }
    }
  /* intensities / raman activities */
  fprintf(f,"[INT]\n");
  for (i=0; i<d->n_vib_modes; i++)
    fprintf(f," %20.10E %20.10E\n",d->freq[i].ir,d->freq[i].raman);
  /* close file */
  file_close(f);
  }

/* -------------------------------------------------------------------------- */
