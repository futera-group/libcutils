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
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* read structure data from the external file
 
   d    - the structure data
   name - name of the file */
void dlpoly_cfg_read(struct dlpoly_cfg *d, char *name) {
  char *line;
  unsigned i;
  FILE *f;
  f = file_open(name,"r");
  /* header */
  d->header = str_read_line_new(f);
  if (!d->header)
    msg_error_f("cannot read header from \"%s\" file",1,name);
  str_trim(d->header);
  /* data ID, pbc ID, number of atoms */
  line = str_read_line_new(f);
  if (!line)
    msg_error_f("cannot read data IDs from \"%s\" file",1,name);
  if (sscanf(line,"%hd%hd%u",&(d->data),&(d->pbc),&(d->n_atoms))!=3)
    msg_error_f("invalid format of data IDs",1);
  line = str_free(line);
  /* cell vectors */
  if (d->pbc) {
    d->cell = mat_falloc(3,3);
    for (i=0; i<3; i++) {
      line = str_read_line_new(f);
      if (!line)
        msg_error_f("cannot read data cell vector #%d",1,i+1);
      if (sscanf(line,"%lf%lf%lf",&(d->cell[i][0]),
          &(d->cell[i][1]),&(d->cell[i][2]))!=3)
        msg_error_f("invalid format of cell vector #%d",1,i+1);
      line = str_free(line);
      }
    }
  /* atomic data */
  d->sym = vec_salloc(d->n_atoms,8);
  d->crd = mat_falloc(d->n_atoms,3);
  if (d->data>0)
    d->vel = mat_falloc(d->n_atoms,3);
  if (d->data>1)
    d->frc = mat_falloc(d->n_atoms,3);
  for (i=0; i<d->n_atoms; i++) {
    line = str_read_line_new(f);
    if (!line)
      msg_error_f("cannot read data atom #%d symbol",1,i+1);
    if (sscanf(line,"%s",d->sym[i])!=1)
      msg_error_f("invalid format of atom #%d symbol",1,i+1);
    line = str_free(line);
    line = str_read_line_new(f);
    if (!line)
      msg_error_f("cannot read data atom #%d coordinates",1,i+1);
    if (sscanf(line,"%lf%lf%lf",&(d->crd[i][0]),
        &(d->crd[i][1]),&(d->crd[i][2]))!=3)
      msg_error_f("invalid format of atom #%d coordinates",1,i+1);
    line = str_free(line);
    if (d->data>0) {
      line = str_read_line_new(f);
      if (!line)
        msg_error_f("cannot read data atom #%d velocity",1,i+1);
      if (sscanf(line,"%lf%lf%lf",&(d->vel[i][0]),
          &(d->vel[i][1]),&(d->vel[i][2]))!=3)
        msg_error_f("invalid format of atom #%d velocity",1,i+1);
      line = str_free(line);
      }
    if (d->data>1) {
      line = str_read_line_new(f);
      if (!line)
        msg_error_f("cannot read data atom #%d force",1,i+1);
      if (sscanf(line,"%lf%lf%lf",&(d->frc[i][0]),
          &(d->frc[i][1]),&(d->frc[i][2]))!=3)
        msg_error_f("invalid format of atom #%d force",1,i+1);
      line = str_free(line);
      }
    }
  file_close(f);
  }

/* -------------------------------------------------------------------------- */
