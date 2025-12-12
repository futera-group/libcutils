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
#include "mol/atom.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* write molecular data to file in xyz format with specific number format
 
   x    - pointer to xyz molecular data struct 
   frmt - format for number printing
   file - pointer to open file */
void mol_xyz_fwritef(struct xyz_mol *x, char *frmt, FILE *file) {
  unsigned i,j;
  if (x->title)
    fprintf(file,"%5d\n%s\n",x->n_atoms,x->title);
  else
    fprintf(file,"%5d\n\n",x->n_atoms);
  for (i=0; i<x->n_atoms; i++) {
    fprintf(file," %-3s",atom_name(x->atom[i].num));
    for (j=0; j<3; j++)
      fprintf(file,frmt,x->atom[i].coord[j]);
    fprintf(file,"\n");
    }
  }

/* write molecular data to file in xyz format 
 
   x    - pointer to xyz molecular data struct 
   file - pointer to open file */
void mol_xyz_fwrite(struct xyz_mol *x, FILE *file) {
  unsigned i,j;
  if (x->title)
    fprintf(file,"%5d\n%s\n",x->n_atoms,x->title);
  else
    fprintf(file,"%5d\n\n",x->n_atoms);
  for (i=0; i<x->n_atoms; i++) {
    fprintf(file," %-3s",atom_name(x->atom[i].num));
    for (j=0; j<3; j++)
      fprintf(file," %18.10f",x->atom[i].coord[j]);
    fprintf(file,"\n");
    }
  }

/* write molecular data to file in xyz format
 
   x    - pointer to xyz molecular data struct 
   name - name of the file */
void mol_xyz_write(struct xyz_mol *x, char *name) {
  FILE *file = stdout;
  if (name && name[0])
    file = file_open(name,"w");
  mol_xyz_fwrite(x,file);
  if (name && name[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
