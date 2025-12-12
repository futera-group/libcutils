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
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Write structure to open GRO file

   d - gromacs data structure
   f - open file stream */
void gmx_dat_fwrite_gro(struct gmx_dat *d, FILE *f) {
  unsigned i,j,k,l,ia,ir;
  /* title */
  if (d->title)
    fprintf(f,"%s\n",d->title);
  else
    fprintf(f,"\n");
  /* number of atoms */
  fprintf(f,"%d\n",d->n_atoms);
  /* structure */
  for (i=0,ia=0,ir=0; i<d->n_frags; i++) {
    for (j=0; j<d->frag[i].n_rep; j++) {
      for (k=0; k<d->frag[i].mol->n_resids; k++,ir++) {
        for (l=0; l<d->frag[i].mol->res[k].n_atoms; l++,ia++) {
          fprintf(f,"%5d%-4s%6s%5d%8.3f%8.3f%8.3f\n",
            ir+1,d->frag[i].mol->res[k].name,
            d->frag[i].mol->res[k].atom[l].name,ia+1,
            d->crd[ia][0]/10.0,d->crd[ia][1]/10.0,d->crd[ia][2]/10.0);
          }
        }
      }
    }
  /* simulation box */
  if (d->box)
    fprintf(f,"%10.5f %10.5f %10.5f\n",
      d->box[0]/10.0,d->box[1]/10.0,d->box[2]/10.0);
  }

/* Write structure to GRO file

   d    - gromacs data structure
   name - name of the file */
void gmx_dat_write_gro(struct gmx_dat *d, char *name) {
  FILE *f = stdout;
  if (name && name[0]) 
    f = file_open(name,"w");
  gmx_dat_fwrite_gro(d,f);
  if (name && name[0]) 
    file_close(f);
  }

/* -------------------------------------------------------------------------- */
