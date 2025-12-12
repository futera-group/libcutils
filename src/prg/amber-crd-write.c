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
#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

/* write one structure to the amber trajectory file

   x    - array with cartesian coordinates
   n    - number of atoms
   box  - box parameters (not writen if NULL)
   file - pointer to open trajectory file */
void amber_crd_write_one(double *x, unsigned n, double *box, FILE *file) {
  unsigned i;
  for (i=0; i<3*n; i++) {
    fprintf(file,"%8.3f",x[i]);
    if ((i+1)%10==0)
      fprintf(file,"\n");
    }
  if ((3*n)%10)
    fprintf(file,"\n");
  if (box) {
    for (i=0; i<3; i++)
      fprintf(file,"%8.3f",box[i]);
    fprintf(file,"\n");
    }
  }

/* -------------------------------------------------------------------------- */

/* write coordinates to amber restart file

   x - vector for cartesian coordinates
   n - number of atoms 
   b - box parameters in case of PBC (not used if NULL)
   b - velocity vector (not used if NULL)
   f - name of the amber restart file */
void amber_crd_write_rst(double *x, unsigned n, double *b, double *v, char *f) {
  unsigned i,j;
  FILE *file = stdout;
  if (f)
    file = file_open(f,"w");
  fprintf(file,"\n%d\n",n);
  /* write coordinates */
  for (i=0; i<n; i++) {
    for (j=0; j<3; j++)
      fprintf(file,"%12.7f",x[3*i+j]);
    if ((i+1)%2==0)
      fprintf(file,"\n");
    }
  if (n%2!=0)
    fprintf(file,"\n");
  /* write velocities */
  if (v) {
    for (i=0; i<n; i++) {
      for (j=0; j<3; j++)
        fprintf(file,"%12.7f",v[3*i+j]);
      if ((i+1)%2==0)
        fprintf(file,"\n");
      }
    if (n%2!=0)
      fprintf(file,"\n");
    }
  /* write box info */
  if (b) {
    for (i=0; i<3; i++)
      fprintf(file,"%12.7f",b[i]);
    for (i=0; i<3; i++)
      fprintf(file,"%12.7f",90.0);
    fprintf(file,"\n");
    }
  if (f)
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
