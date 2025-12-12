/******************************************************************************\
 *                                                                            * 
 *  Libcutils - library of C function                                         * 
 *                                                                            *
 *  Version:             3.4                                                  * 
 *  Date:                26/03/2020                                           *
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
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* file format */
#define FRMT_FULL   1  /* full data set (res-name,atom-name,id,coords) */
#define FRMT_COORD  2  /* coordinates only */

/* -------------------------------------------------------------------------- */

/* find specified data block in G96 file

   f   - open file stream
   ter - terminate program in case of error or not
   key - keyword to search for */
int gmx_dat_fread_g96_block(FILE *f, short ter, char *key) {
  char *line;
  line = str_ffind_b_new(f,key);
  if (!line) {
    if (!ter)
      return(0);
    msg_error_f("cannot find '%s' data block in G96 file",1,key);
    }
  str_free(line);
  return(1);
  }

/* Read structure from open G96 file
  
   d     - gromacs data structure
   step  - simulation step (output)
   time  - simulation time (output)
   title - check whether there is a title or not
   ter   - terminate program in case of error or not
   f     - open file stream */
int gmx_dat_fread_g96(struct gmx_dat *d, long unsigned *step, double *time, 
  short *title, short ter, FILE *f) {
  unsigned i,j,k,l,n,ia;
  long unsigned st;
  double tm,crd[3];
  char *line,**w;
  short frmt = 0;
  fpos_t fpos;
  /* title */
  fgetpos(f,&fpos);
  line = str_read_line_new(f);
  if (strstr(line,"TITLE")) {
    (*title) = 1;
    d->title = str_read_line_new(f);
    str_trim(d->title);
    }
  else {
    (*title) = 0;
    fsetpos(f,&fpos);
    }
  str_free(line);
  /* time step (optional) */
  if (gmx_dat_fread_g96_block(f,0,"TIMESTEP")) {
    line = str_read_line_new(f);
    if (sscanf(line,"%lu%lf",&st,&tm)!=2) {
      if (!ter)
        return(0);
      msg_error("cannot read time step in G96 file",1);
      }
    if (step)
      (*step) = st;
    if (time)
      (*time) = tm;
    str_free(line);
    }
  else
    rewind(f);
  /* atoms */
  if (!gmx_dat_fread_g96_block(f,ter,"POSITION"))
    return(0);
  d->crd = mat_ffree(d->crd,d->n_atoms);
  d->crd = mat_falloc(d->n_atoms,3);
  for (i=0,ia=0; i<d->n_frags; i++)
    for (j=0; j<d->frag[i].n_rep; j++)
      for (k=0; k<d->frag[i].mol->n_resids; k++) {
        for (l=0; l<d->frag[i].mol->res[k].n_atoms; l++,ia++) {
          line = str_read_line_new(f);
          if (!line) {
            if (!ter)
              return(0);
            msg_error("unexpected end of G96 file"
              " while reading coordinates",1);
            }
          /* check format */
          if (ia==0) {
            w = str_split(line,' ',&n);
            if (n >= 7)
              frmt = FRMT_FULL;
            else
              frmt = FRMT_COORD;
            vec_sfree(w,n);
            }
          /* read data */
          switch (frmt) {
            case FRMT_FULL:
              if (sscanf(line,"%*d%*s%*s%*d%lf%lf%lf",
                    &(crd[0]),&(crd[1]),&(crd[2]))!=3) {
                if (!ter)
                  return(0);
                msg_error_f("cannot read coordinates of"
                  " atom #%d in G96 file",1,i+1);
                }
              break;
            case FRMT_COORD:
              if (sscanf(line,"%lf%lf%lf",&(crd[0]),&(crd[1]),&(crd[2]))!=3) {
                if (!ter)
                  return(0);
                msg_error_f("cannot read coordinates of"
                  " atom #%d in G96 file",1,i+1);
                }
              break;
            }
          vec_fcopy_scaled(d->crd[ia],crd,10.0,3);
          str_free(line);
          }
        }
  /* simulation box */
  if (!gmx_dat_fread_g96_block(f,ter,"BOX"))
    return(0);
  line = str_read_line_new(f);
  if (line && sscanf(line,"%lf%lf%lf",&(crd[0]),&(crd[1]),&(crd[2]))==3) {
    d->box = vec_ffree(d->box);
    d->box = vec_fcopy_new(crd,3);
    vec_fscale(d->box,10.0,3);
    }
  return(1);
  }

/* Read structure from G96 file

   d    - gromacs data structure
   name - name of the file */
void gmx_dat_read_g96(struct gmx_dat *d, char *name) {
  short title;
  FILE *f;
  f = file_open(name,"r");
  gmx_dat_fread_g96(d,NULL,NULL,&title,1,f);
  file_close(f);
  }

/* -------------------------------------------------------------------------- */
