/******************************************************************************\
 *                                                                            * 
 *  Libcutils - library of C function                                         * 
 *                                                                            *
 *  Version:             3.4                                                  * 
 *  Date:                07/11/2021                                           *
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
#include <cmn/message.h>
#include <cmn/string.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* read header of the Lammps DCD trajectory file

   d    - system data 
   name - name of the DCD file */
FILE* lammps_trj_dcd_read_header(struct lammps_dat *d, char *name) {
  unsigned i,j,k,n_titles;
  char s[80];
  int n[41];
  FILE *f;
  union {int n; float r; char s[4];} u;
  /* open the file */
  f = file_open(name,"rb");
  /* leading array of config IDs */
  file_bin_frt_read_int(n,21,f);
  /* info flag */
  u.n = n[0];
  for (i=0; i<4; i++)
    s[i] = u.s[i];
  s[i] = '\0';
  if (!str_compare(s,"CORD"))
    msg_error_f("invalid format of DCD file (\"%s\" label)",1,s);
  /* memory allocation */
  if (!d->trj)
    d->trj = lammps_trj_new();
  /* number of frames */
  d->trj->n_frames = n[1];
  d->trj->save_freq = n[3];
  d->trj->n_steps = n[4];
  /* timestep */
  u.n = n[10];
  d->trj->timestep = u.r;
  /* simulation box */
  d->trj->has_box = n[11];
  /* dimension */
  d->trj->has_4d = n[12];
  if (d->trj->has_4d)
    msg_error("4D data found in the DCD file",1);
  /* fluctuation charges */
  d->trj->has_chrg = n[13];
  /* title & remarks */
  file_bin_frt_read_int(n,41,f);
  n_titles = n[0];
  if (n_titles < 1)
    msg_error("no title line found in DCD file",1);
  if (n_titles > 2)
    msg_error_f("unexpected number of title lines (%d)",1,n_titles);
  for (i=0; i<n_titles; i++) {
    for (j=0; j<20; j++) {
      u.n = n[20*i+j+1];
      for (k=0; k<4; k++)
        d->trj->title[i][4*j+k] = u.s[k];
      }
    d->trj->title[i][80] = '\0';
    }
  /* number of atoms */
  file_bin_frt_read_int(n,1,f);
  if (n[0] != d->n_atoms)
    msg_error_f("inconsistent number of atoms in DCD file (%d/%d)",1,
      n[0],d->n_atoms);
  return(f);
  }

/* read one structure from Lammps DCD trajectory file
 
   d - system data
   v - array for reading N single-precision coordinates
   f - open file stream */
int lammps_trj_dcd_read_frame(struct lammps_dat *d, float *v, FILE *f) {
  unsigned i;
  double r[6];
  /* cell size */
  if (d->trj->has_box) {
    file_bin_frt_read_real(r,6,f);
    lammps_box_set(d->box,r[0],r[2],r[5],r[4],r[3],r[1]);
    }
  /* coordinates */
  file_bin_frt_read_real_s(v,d->n_atoms,f);
  for (i=0; i<d->n_atoms; i++)
    d->atom[i].crd[0] = v[i];
  file_bin_frt_read_real_s(v,d->n_atoms,f);
  for (i=0; i<d->n_atoms; i++)
    d->atom[i].crd[1] = v[i];
  file_bin_frt_read_real_s(v,d->n_atoms,f);
  for (i=0; i<d->n_atoms; i++)
    d->atom[i].crd[2] = v[i];
  return(1);
  }

/* -------------------------------------------------------------------------- */
