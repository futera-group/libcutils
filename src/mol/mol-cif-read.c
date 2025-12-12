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
#include <cmn/message.h>
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/types.h>
#include "mol/atom.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* read one keyword from cif file
 
   w - name of the keyword
   t - type of data value
   d - data storage for keyword value
   f - pointer to open cif file */
void mol_cif_fread_flag(char *w, short t, void *d, FILE *f) {
  char *line,val[256];
  unsigned id;
  line = str_ffind_new(f,w);
  if (!line)
    msg_error_f("keyword \"%s\" not found in CIF file",1,w);
  /* string value */
  if (t==TYPE_STRING) {
    str_trim(line);
    if (!str_char_first_b(line,&id))
      msg_error_f("invalid format of keyword \"%s\"",1,w);
    (*((char**)d)) = str_copy_new(line+id);
    str_trim(*((char**)d));
    }
  /* numberic values */
  else {
    if (sscanf(line,"%*s%s",val)!=1)
      msg_error_f("invalid format of keyword \"%s\"",1,w);
    if (!type_read(val,d,t))
      msg_error_f("cannot read value of keyword \"%s\"",1,w);
    }
  str_free(line);
  }

/* read molecular data from already open cif file
 
   c - pointer to cif molecular data struct
   f - pointer to open cif file */
void mol_cif_fread(struct cif_mol *c, FILE *f) {
  char *line,atom_name[80],s1[80],s2[80];
  double x,cell_dim[3],cell_angle[3];
  unsigned id = 0;
  struct cif_atom *a;
  struct queue *q;
  /* cell data */
  mol_cif_fread_flag("_cell_length_a",TYPE_DOUBLE,
    &(cell_dim[0]),f);
  mol_cif_fread_flag("_cell_length_b",TYPE_DOUBLE,
    &(cell_dim[1]),f);
  mol_cif_fread_flag("_cell_length_c",TYPE_DOUBLE,
    &(cell_dim[2]),f);
  mol_cif_fread_flag("_cell_angle_alpha",TYPE_DOUBLE,
    &(cell_angle[0]),f);
  mol_cif_fread_flag("_cell_angle_beta",TYPE_DOUBLE,
    &(cell_angle[1]),f);
  mol_cif_fread_flag("_cell_angle_gamma",TYPE_DOUBLE,
    &(cell_angle[2]),f);
  cell_set_side_angle_v(c->cell,cell_dim,cell_angle);
  /* space group */
  mol_cif_fread_flag("_symmetry_space_group_name_H-M",TYPE_STRING,
    &(c->cell->space_group_name),f);
  mol_cif_fread_flag("_symmetry_int_tables_number",TYPE_UINT,
    &(c->cell->space_group_id),f);
  /* atoms */
  line = str_ffind_new(f,"_atom_site_type_symbol");
  if (!line)
    msg_error("atom specificiation missing in CIF file",1);
  str_free(line);
  /* read data to queue */
  q = queue_alloc();
  for (line=str_read_line_new(f); line;
       str_free(line),line=str_read_line_new(f)) {
    a = mol_cif_atom_new(1);
    if (sscanf(line,"%s%lf%lf%lf%lf%s%lf%s",atom_name,&x,
      &(a->coord[0]),&(a->coord[1]),&(a->coord[2]),s1,&x,s2)!=8)
      break;
    a->name = str_copy_new(atom_name);
    a->num = atom_num(s2);
    queue_add(q,a);
    }
  /* convert queue to array */
  c->n_atoms = q->num;
  c->atom = mol_cif_atom_new(c->n_atoms);
  while (q->num) {
    a = queue_get(q);
    mol_cif_atom_copy(a,c->atom+id);
    mol_cif_atom_free(a,1);
    id++;
    }
  /* clean memory */
  queue_free(q);
  }

/* read molecular data from cif file
 
   c - pointer to cif molecular data struct
   f - name of the file */
void mol_cif_read(struct cif_mol *c, char *f) {
  FILE *file = file_open(f,"r");
  mol_cif_fread(c,file);
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
