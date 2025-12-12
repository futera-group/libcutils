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
#include <stdlib.h>
#include <cmn/file.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* read molecular data from open cube file
 
   c - pointer to cube molecular data struct
   f - pointer to open file */
void mol_cub_fread(struct cub_mol *c, FILE *f) {
  char *line;
  unsigned i,j,k;
  int n;
  /* memory allocation */
  if (!c->grid)
    c->grid = mol_cub_grid_new();
  /* title */
  c->title = str_read_line_new(f);
  if (!c->title)
    msg_error("cube file is empty",1);
  str_rtrim(c->title);
  /* comment line */
  c->desc = str_read_line_new(f);
  if (!c->desc)
    msg_error("unexpected error while reading header of cube file",1);
  str_rtrim(c->desc);
  /* number of atoms and grid origin */
  line = str_read_line_new(f);
  if (!line)
    msg_error("unexpected error while reading origin of mesh",1);
  if (sscanf(line,"%d%lf%lf%lf",&n,&(c->grid->origin[0]),
    &(c->grid->origin[1]),&(c->grid->origin[2]))!=4)
    msg_error("invalid format of mesh origin",1);
  c->n_atoms = abs(n);
  str_free(line);
  /* grid vectors */
  for (i=0; i<3; i++) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected error while reading grid vectors",1);
    if (sscanf(line,"%u%lf%lf%lf",&(c->grid->range[i]),&(c->grid->vector[i][0]),
      &(c->grid->vector[i][1]),&(c->grid->vector[i][2]))!=4)
      msg_error("invalid format of grid vectors",1);
    str_free(line);
    }
  /* atomic coordinates */
  c->atom = mol_cub_atom_new(c->n_atoms);
  for (i=0; i<c->n_atoms; i++) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected error while reading atom specification",1);
    if (sscanf(line,"%u%*f%lf%lf%lf",&(c->atom[i].num),&(c->atom[i].coord[0]),
      &(c->atom[i].coord[1]),&(c->atom[i].coord[2]))!=4)
      msg_error("invalid format of atom specification",1);
    str_free(line);
    }
  /* grid data */
  if (!c->data)
    c->data = mol_cub_data_new(c->grid);
  for (i=0; i<c->grid->range[0]; i++)
    for (j=0; j<c->grid->range[1]; j++)
      for (k=0; k<c->grid->range[2]; k++)
        if (fscanf(f,"%lf",&(c->data[i][j][k]))!=1)
          msg_error("unexpected error while reading grid data",1);
  }

/* read molecular data from cube file
 
   c - pointer to cube molecular data struct
   f - name of the file */
void mol_cub_read(struct cub_mol *c, char*f) {
  FILE *file;
  file = file_open(f,"r");
  mol_cub_fread(c,file);
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
