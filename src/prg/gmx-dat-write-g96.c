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

/* Write structure to open G96 file

   d     - gromacs data structure
   step  - simulation step (output)
   time  - simulation time (output)
   title - write title or not
   f     - open file stream */
void gmx_dat_fwrite_g96(struct gmx_dat *d, long unsigned step, double time,
  short title, FILE *f) {
  unsigned i,j,k,l,n,ia;
  /* title */
  if (title) {
    fprintf(f,"TITLE\n");
    if (d->title)
      fprintf(f,"%s\n",d->title);
    fprintf(f,"END\n");
    }
  /* time step */
  fprintf(f,"TIMESTEP\n");
  fprintf(f,"%15ld%15.6f\n",step,time);
  fprintf(f,"END\n");
  /* coordinates */
  fprintf(f,"POSITIONRED\n");
  for (i=0,ia=0; i<d->n_frags; i++)
    for (j=0; j<d->frag[i].n_rep; j++)
      for (k=0; k<d->frag[i].mol->n_resids; k++) {
        for (l=0; l<d->frag[i].mol->res[k].n_atoms; l++,ia++) {
          for (n=0; n<3; n++)
            fprintf(f,"%15.9f",d->crd[ia][n]/10.0);
          fprintf(f,"\n");
          }
        }
  fprintf(f,"END\n");
  /* simulation box */
  if (d->box) {
    fprintf(f,"BOX\n");
    for (i=0; i<3; i++)
      fprintf(f,"%15.9f",d->box[i]/10.0);
    fprintf(f,"\n");
    fprintf(f,"END\n");
    }
  }

/* Write structure to G96 file

   d    - gromacs data structure
   name - name of the file */
void gmx_dat_write_g96(struct gmx_dat *d, char *name) {
  FILE *f = stdout;
  if (name && name[0]) 
    f = file_open(name,"w");
  gmx_dat_fwrite_g96(d,0,0.0,1,f);
  if (name && name[0]) 
    file_close(f);
  }

/* -------------------------------------------------------------------------- */
