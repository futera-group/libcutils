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
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/units.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include <mol/molec.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* read one line from gaussian input file
 
   ne - storage for line data
   what - short description of expected data
   file - pointer to open file */
char* gauss_inp_read_line(char *what, FILE *file) {
  char *line;
  line = str_read_line_new(file);
  if (!line)
    msg_error_f("unexpected end of gaussian input file while reading %s",
      1,what);
  return(line);
  }

/* read background charges from open gaussian input file
 
   g - pointer to gaussian data struct
   f - name of the input file */
void gauss_inp_read_chrg_f(struct gauss_dat *g, FILE *f) {
  char *line;
  unsigned i,j,id;
  double *h;
  struct queue *q;
  /* update charges */
  if (g->n_bg_charges && g->bg_chrg) {
    for (i=0; i<g->n_bg_charges; i++) {
      line = gauss_inp_read_line("background charges",f);
      if (sscanf(line,"%lf%lf%lf%lf",&(g->bg_chrg[i][0]),
          &(g->bg_chrg[i][1]),&(g->bg_chrg[i][2]),&(g->bg_chrg[i][3]))!=4)
        msg_error_f("invalid format of background charges in %s file",1,f);
      for (j=0; j<3; j++)
        g->bg_chrg[i][j] /= CONV_B_ANG;
      line = str_free(line);
      }
    }
  /* new array of charges */
  else {
    /* read data to queue */
    q = queue_alloc();
    for (line=str_read_line_new(f); line;
         line=str_free(line),line=str_read_line_new(f)) {
      str_trim(line);
      if (line[0]=='\0' || line[0]=='\n')
        break;
      h = vec_falloc(4);
      if (sscanf(line,"%lf%lf%lf%lf",&(h[0]),&(h[1]),&(h[2]),&(h[3]))!=4)
        msg_error("invalid format of charges in gaussian input file",1);
      queue_add(q,h);
      }
    /* convert queue to gaussian format */
    g->n_bg_charges = q->num;
    g->bg_chrg = mat_falloc(g->n_bg_charges,4);
    id = 0;
    while (q->num) {
      h = queue_get(q);
      for (i=0; i<3; i++)
        g->bg_chrg[id][i] = h[i]/CONV_B_ANG;
      g->bg_chrg[id][3] = h[3];
      id++;
      }
    queue_free(q);
    }
  }

/* read background charges from gaussian input file
 
   g - pointer to gaussian data struct
   f - name of the input file */
void gauss_inp_read_chrg(struct gauss_dat *g, char *f) {
  char *line;
  unsigned i;
  FILE *file;
  /* open file */
  file = file_open(f,"r");
  /* block of instructions */
  for (line=str_read_line_new(file); line;
       line=str_free(line),line=str_read_line_new(file)) {
    str_trim(line);
    if (line[0]=='\0' || line[0]=='\n')
      break;
    }
  /* comment, change and multiplicity */
  for (i=0; i<3; i++) {
    line = gauss_inp_read_line("comment",file);
    line = str_free(line);
    }
  /* coordinates */
  for (line=str_read_line_new(file); line;
       line=str_free(line),line=str_read_line_new(file)) {
    str_trim(line);
    if (line[0]=='\0' || line[0]=='\n')
      break;
    }
  /* point charges */
  gauss_inp_read_chrg_f(g,file);
  /* close file */
  file_close(file);
  }

/* read coordinates from open gaussian input file
 
   g - pointer to gaussian data struct
   f - pointer to open file stream */
void gauss_inp_read_coord_f(struct gauss_dat *g, FILE *f) {
  char *line,sym[80];
  unsigned i,j,id;
  struct xyz_atom *a;
  struct queue *q;
  /* update coordinates */
  if (g->n_atoms && g->atom) {
    for (i=0; i<g->n_atoms; i++) {
      line = gauss_inp_read_line("cooordinates",f);
      if (sscanf(line,"%*s%lf%lf%lf",&(g->atom[i].coord[0]),
          &(g->atom[i].coord[1]),&(g->atom[i].coord[2]))!=3)
        msg_error_f("invalid format of coordinates in %s file",1,f);
      for (j=0; j<3; j++)
        g->atom[i].coord[j] /= CONV_B_ANG;
      line = str_free(line);
      }
    }
  /* new coordinates */
  else {
    /* read data to queue */
    q = queue_alloc();
    for (line=str_read_line_new(f); line;
         line=str_free(line),line=str_read_line_new(f)) {
      str_trim(line);
      if (line[0]=='\0' || line[0]=='\n')
        break;
      a = mol_xyz_atom_new(1);
      if (sscanf(line,"%s%lf%lf%lf",sym,&(a->coord[0]),&(a->coord[1]),
          &(a->coord[2]))!=4)
        msg_error("invalid format of coordinates in gaussian input file",1);
      a->num = atom_num(sym);
      queue_add(q,a);
      }
    /* convert queue to gaussian format */
    g->n_atoms = q->num;
    g->atom = gauss_atom_new(g->n_atoms);
    id = 0;
    while (q->num) {
      a = queue_get(q);
      g->atom[id].num = a->num;
      g->atom[id].charge[GAUSS_CH_NC] = (double)g->atom[id].num;
      for (i=0; i<3; i++)
        g->atom[id].coord[i] = (a->coord[i]/=CONV_B_ANG);
      id++;
      }
    queue_free(q);
    }
  }

/* read coordinates from gaussian input file
 
   g - pointer to gaussian data struct
   f - name of the input file */
void gauss_inp_read_coord(struct gauss_dat *g, char *f) {
  char *line;
  unsigned i;
  FILE *file;
  /* open file */
  file = file_open(f,"r");
  /* block of instructions */
  for (line=str_read_line_new(file); line;
       line=str_free(line),line=str_read_line_new(file)) {
    str_trim(line);
    if (line[0]=='\0' || line[0]=='\n')
      break;
    }
  /* comment, change and multiplicity */
  for (i=0; i<3; i++) {
    line = gauss_inp_read_line("comment",file);
    line = str_free(line);
    }
  /* coordinates */
  gauss_inp_read_coord_f(g,file);
  /* close file */
  file_close(file);
  }

/* read instructions, charge and multiplicity from open gaussian input file
 
   g - pointer to gaussian data struct
   f - pointer to open file stream */
void gauss_inp_read_header_f(struct gauss_dat *g, FILE *f) {
  char *line;
  /* block of instructions */
  for (line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    str_trim(line);
    if (line[0]=='\0' || line[0]=='\n')
      break;
    }
  /* comment */
  g->job_title = gauss_inp_read_line("comment",f);
  str_trim(g->job_title);
  /* change and multiplicity */
  line = gauss_inp_read_line("second free line",f);
  line = str_free(line);
  line = gauss_inp_read_line("charge and multiplicity",f);
  if (sscanf(line,"%d%u",&(g->charge),&(g->multiplicity))!=2)
    msg_error("invalid format of charge and multiplicity"
     " in gaussian input file",1);
  line = str_free(line);
  }

/* read instructions, charge and multiplicity from gaussian input file
 
   g - pointer to gaussian data struct
   f - name of the input file */
void gauss_inp_read_header(struct gauss_dat *g, char *f) {
  FILE *file;
  file = file_open(f,"r");
  gauss_inp_read_header_f(g,file);
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
