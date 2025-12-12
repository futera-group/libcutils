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
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* write structure data to external file
 
   d    - the structure data
   name - name of the file */
void dlpoly_cfg_write(struct dlpoly_cfg *d, char *name) {
  unsigned i,j;
  FILE *f = stdout;
  if (name && name[0]) 
    f = file_open(name,"w");
  /* header */
  fprintf(f,"%s\n",d->header);
  /* data ID, pbc ID, number of atoms */
  fprintf(f,"%10d%10d%10d\n",d->data,d->pbc,d->n_atoms);
  /* cell vectors */
  if (d->pbc) {
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++)
        fprintf(f,"%20.12f",d->cell[i][j]);
      fprintf(f,"\n");
      }
    }
  /* atomic data */
  for (i=0; i<d->n_atoms; i++) {
    fprintf(f,"%-8s%10d\n",d->sym[i],i+1);
    fprintf(f,"%20.10f%20.10f%20.10f\n",
      d->crd[i][0],d->crd[i][1],d->crd[i][2]);
    if (d->data>0)
      fprintf(f,"%20.10f%20.10f%20.10f\n",
        d->vel[i][0],d->vel[i][1],d->vel[i][2]);
    if (d->data>1)
      fprintf(f,"%20.10f%20.10f%20.10f\n",
        d->frc[i][0],d->frc[i][1],d->frc[i][2]);
    }
  if (name && name[0]) 
    file_close(f);
  }

/* -------------------------------------------------------------------------- */
