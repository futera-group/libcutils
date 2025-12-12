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
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* write grid parameters to file in cube format
 
   c - pointer to cube molecular data struct
   n - number of atoms
   f - pointer to open file */
void mol_cub_fwrite_parm(struct cub_grid *g, unsigned n, FILE *f) {
  unsigned i;
  fprintf(f,"%5d%12.6f%12.6f%12.6f\n",n,g->origin[0],g->origin[1],g->origin[2]);
  for (i=0; i<3; i++)
    fprintf(f,"%5d%12.6f%12.6f%12.6f\n",g->range[i],g->vector[i][0],
      g->vector[i][1],g->vector[i][2]);
  }

/* write molecular data to file in cube format
 
   c - pointer to cube molecular data struct
   f - pointer to open file */
void mol_cub_fwrite(struct cub_mol *c, FILE *f) {
  unsigned i,j,k;
  fprintf(f,"%s\n",c->title);
  fprintf(f,"%s\n",c->desc);
  mol_cub_fwrite_parm(c->grid,c->n_atoms,f);
  for (i=0; i<c->n_atoms; i++)
    fprintf(f,"%5d%12.6f%12.6f%12.6f%12.6f\n",c->atom[i].num,
      (double)(c->atom[i].num),c->atom[i].coord[0],
      c->atom[i].coord[1],c->atom[i].coord[2]);
  for (i=0; i<c->grid->range[0]; i++)
    for (j=0; j<c->grid->range[1]; j++) {
      for (k=0; k<c->grid->range[2]; k++) {
        fprintf(f,"%13.5E",c->data[i][j][k]);
        if (!((k+1)%6))
          fprintf(f,"\n");
        }
      if (c->grid->range[2]%6)
        fprintf(f,"\n");
      }
  }

/* write molecular data to file in cube format
 
   c - pointer to cube molecular data struct
   f - name of the file */
void mol_cub_write(struct cub_mol *c, char *f) {
  FILE *file = stdout;
  if (f && f[0])
    file = file_open(f,"w");
  mol_cub_fwrite(c,file);
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
