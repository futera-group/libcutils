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
#include <string.h>
#include <cmn/file.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

#define AMBER_CRD_NCOL 10  /* number of values in one line */
#define AMBER_CRD_NUMW  8  /* width of one number in the line */
#define AMBER_CRD_NUMD  4  /* position of decimal point in a number */

/* -------------------------------------------------------------------------- */

/* read numbers from one line in amber trajectory file (non-standard, overfull)
 
   x    - array of numbers
   n    - ID of last number in the array
   file - pointer to open trajectory file */
int amber_crd_read_line_nstd(double *x, unsigned *n, unsigned *l, FILE *file) {
  char *line,snum[256];
  unsigned i,m,w,num = 0,off = 0,id=(*n);
  /* read line */
  line = str_read_line_new(file);
  if (!line)
    return(num);
  (*l) = (*l)+1;
  /* read numbers */
  m = str_length(line);
  for (num=0; num<AMBER_CRD_NCOL; num++) {
    strncpy(snum,line+off+AMBER_CRD_NUMW*num,AMBER_CRD_NUMW);
    snum[AMBER_CRD_NUMW]='\0';
    /* invalid number format */
    if (snum[AMBER_CRD_NUMD]!='.') {
      i = off+AMBER_CRD_NUMW*num;
      while (i<m && line[i]!='.') i++;
      if (i<m) {
        w = i-off-AMBER_CRD_NUMW*num+4;
        strncpy(snum,line+off+AMBER_CRD_NUMW*num,w);
        off += (w-AMBER_CRD_NUMW);
        snum[w]='\0';
        }
      }
    if (sscanf(snum,"%lf",&(x[id++]))!=1) {
      return(num);
      }
    }
  str_free(line);
  (*n) += num;
  return(num);
  }

/* read numbers from one line in amber trajectory file
 
   x    - array of numbers
   n    - ID of last number in the array
   file - pointer to open trajectory file */
int amber_crd_read_line(double *x, unsigned *n, unsigned *l, FILE *file) {
  char *line,snum[AMBER_CRD_NUMW+1];
  unsigned num = 0,id = (*n);
  /* read line */
  line = str_read_line_new(file);
  if (!line)
    return(num);
  /* read numbers */
  for (num=0; num<AMBER_CRD_NCOL; num++) {
    strncpy(snum,line+AMBER_CRD_NUMW*num,AMBER_CRD_NUMW);
    snum[AMBER_CRD_NUMW]='\0';
    if (sscanf(snum,"%lf",&(x[id++]))!=1) {
      return(num);
      }
    }
  str_free(line);
  (*n) += num;
  return(num);
  }

/* read one structure from amber trajectory file

   x    - array for read coordinates
   n    - number of atoms
   b    - array for box size (no box info read if NULL)
   file - pointer to open trajectory file */
int amber_crd_read_one(double *x, unsigned n, double *b, FILE *file) {
  short rval = 1,nostd = 1;
  unsigned i,r,id = 0,lnum = 1;
  char *line;
  int (*fce_line_read)(double*,unsigned*,unsigned*,FILE*);
  /* function for line reading */
  fce_line_read = (nostd ? amber_crd_read_line_nstd : amber_crd_read_line);
  /* read full lines */
  for (i=0; i<(3*n)/AMBER_CRD_NCOL; i++) 
    if (fce_line_read(x,&id,&lnum,file)!=AMBER_CRD_NCOL) {
      rval = 0;
      break;
      }
  if (rval) {
    /* last incomplete line */
    r = (3*n)%AMBER_CRD_NCOL;
    if (r && fce_line_read(x,&id,&lnum,file)!=r) 
      rval = 0;
    /* box info */
    if (rval && b) {
      line = str_read_line_new(file);
      if (!line)
        rval = 0;
      if (rval && sscanf(line,"%lf%lf%lf",&(b[0]),&(b[1]),&(b[2]))!=3) 
        rval = 0;
      str_free(line);
      }
    }
  return(rval);
  }

/* read specified structure from amber trajectory file

   x    - array for read coordinates
   n    - number of atoms
   b    - array for box size (no box info read if NULL)
   id   - specify the structure for reading
   file - pointer to open trajectory file */
int amber_crd_read_id(double *x, unsigned n, double *b, unsigned id,
  FILE *file) {
  int ret = 1;
  char *line;
  unsigned i;
  /* set file pointer to the beginning */
  rewind(file);
  /* check if trajectory file is not empty */
  line = str_read_line_new(file);
  if (!line)
    ret=0;
  else {
    /* drop first n-1 structures */
    for (i=1; i<(id+1); i++)
      if (!amber_crd_read_one(x,n,b,file)) {
        ret = 0;
        break;
        }
    /* read the specified structure */
    if (ret && !amber_crd_read_one(x,n,b,file)) 
      ret = 0;
    }
  return(ret);
  }

/* -------------------------------------------------------------------------- */

/* read coordinates from amber restart file

   x - vector for cartesian coordinates
   n - number of atoms 
   b - box parameters in case of PBC (not read if NULL)
   v - vector for velocities (not read if NULL)
   f - name of the amber restart file */
short amber_crd_read_rst(double *x, unsigned n, double *b, double *v, char *f) {
  unsigned i,j,nf;
  char *line;
  double *t;
  fpos_t f_pos;
  short ret = 0,stop = 0;
  FILE *file;
  file = file_open(f,"r");
  /* read header */
  line = str_read_line_new(file);
  if (!line)
    msg_error("amber coordinate file is empty",1);
  line = str_free(line);
  line = str_read_line_new(file);
  if (!line)
    msg_error("unexpected end of amber coordinate file",1);
  if (sscanf(line,"%u",&nf)!=1)
    msg_error("cannot read number of atoms in amber coordinate file",1);
  if (nf!=n)
   msg_error("inconsistent number of atoms in amber coordinate file",1);
  line = str_free(line);
  /* read coordinates */
  for (i=0; i<n; i++)
    for (j=0; j<3; j++)
      if (fscanf(file,"%lf",&(x[3*i+j]))!=1)
        msg_error("cannot read coodinates from amber coordinate file",1);
  /* read velocities */
  t = (v ? v : vec_falloc(3*n));
  fgetpos(file,&f_pos);
  for (i=0; !stop && i<n; i++)
    for (j=0; j<3; j++)
      if (fscanf(file,"%lf",&(t[3*i+j]))!=1) {
        fsetpos(file,&f_pos);
        stop = 1;
        break;
        }
  if (v)
    ret = (stop ? 0 : 1);
  else
    t = vec_ffree(t);
  /* read box info */
  if (b) {
    for (i=0; i<3; i++)
      if (fscanf(file,"%lf",&(b[i]))!=1)
        msg_error("cannot read box parameters from amber coordinate file",1);
    }
  file_close(file);
  return(ret);
  }

/* -------------------------------------------------------------------------- */
