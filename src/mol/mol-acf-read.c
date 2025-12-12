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
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "mol/atom.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* read atomic specification from open acf file

   c - pointer to acf molecular data struct
   f - pointer to open acf file */
void mol_acf_fread_atoms(struct acf_mol *c, FILE *f) {
  char *line,res_name[256],atom_name[256],atom_type[256];
  unsigned i,id = 0;
  struct queue *q;
  struct acf_atom *at;
  /* read data from file */
  q = queue_alloc();
  for (line=str_read_line_new(f); line; 
       line=str_free(line),line=str_read_line_new(f)) {
    if (!str_sub_bfind(line,"ATOM"))
      break;
    at = mol_acf_atom_new(1);
    /* atom name */
    for (i=0; i<6; i++)
      atom_name[i] = line[11+i];
    atom_name[6] = '\0';
    str_trim(atom_name);
    /* residuum name */
    for (i=0; i<3; i++)
      res_name[i] = line[17+i];
    res_name[3] = '\0';
    str_trim(res_name);
    /* coordinates */
    if (sscanf(line+26,"%lf%lf%lf",
      &(at->coord[0]),&(at->coord[1]),&(at->coord[2]))!=3)
      msg_error("cannot read coordinate from AC file",1);
    /* charge */
    if (sscanf(line+54,"%lf",&(at->charge))!=1) 
      msg_error("cannot read charge from AC file",1);
    /* atomic type */
    if (sscanf(line+64,"%s",atom_type)!=1) 
      msg_error("cannot read atomic types from AC file",1);
    /* save data */
    at->name = str_copy_new(atom_name);
    at->type = str_copy_new(atom_type);
    at->num = atom_num_pdb(atom_name);
    if (!c->name) 
      c->name = str_copy_new(res_name);
    queue_add(q,at);
    }
  /* convert data to acf structure */
  c->n_atoms = q->num;
  c->atom = mol_acf_atom_new(c->n_atoms);
  c->charge = 0.0;
  while (q->num) {
    at = queue_get(q);
    mol_acf_atom_copy(at,c->atom+id);
    c->charge += at->charge;
    mol_acf_atom_free(at,1);
    id++;
    }
  /* clean memory */
  queue_free(q);
  }

/* read bond specification from open acf file

   c - pointer to acf molecular data struct
   f - pointer to open acf file */
void mol_acf_fread_bonds(struct acf_mol *c, FILE *f) {
  char *line;
  unsigned id = 0,*b;
  struct queue *q;
  rewind(f);
  q = queue_alloc();
  for (line=str_read_line_new(f); line; 
       line=str_free(line),line=str_read_line_new(f)) 
    if (str_sub_bfind(line,"BOND")) {
      b = vec_ualloc(2);
      if (sscanf(line,"%*s%*d%d%d%*d%*s%*s",&(b[0]),&(b[1]))!=2)
        msg_error("invalid format of bond specification",1);
      queue_add(q,b);
      }
  c->n_bonds = q->num;
  c->bond = mat_ualloc(c->n_bonds,2);
  while (q->num) {
    b = queue_get(q);
    c->bond[id][0] = b[0]-1;
    c->bond[id][1] = b[1]-1;
    vec_ufree(b);
    id++;
    }
  /* clean memory */
  queue_free(q);
  }

/* read molecular data from open acf file

   c - pointer to acf molecular data struct
   f - pointer to open acf file */
void mol_acf_fread(struct acf_mol *c, FILE *f) {
  unsigned id;
  char *line;
  /* charge */
  line = str_read_line_new(f);
  if (!line)
    msg_error("acf molecular file is empty",1);
  if (!str_sub_bfind(line,"CHARGE"))
    msg_error("charge keyword expected in the first line",1);
  if (sscanf(line,"%*s %lf",&(c->charge))!=1)
    msg_error("invalid format of charge specification",1);
  str_free(line);
  /* formula */
  line = str_read_line_new(f);
  if (!line)
    msg_error("unexpected end of acf file while reading formula",1);
  if (!str_sub_bfind(line,"Formula:"))
    msg_error("formula specification expected in the second line",1);
  str_trim(line);
  if (!str_char_first_b(line,&id))
    msg_error_f("invalid formula specification",1);
  c->form = str_copy_new(line+id);
  str_trim(c->form);
  str_free(line);
  /* atomic coordinates */
  mol_acf_fread_atoms(c,f);
  /* interatomic bonds */
  mol_acf_fread_bonds(c,f);
  }

/* read molecular data from acf file
 
   c - pointer to cif molecular data struct
   f - name of the file */
void mol_acf_read(struct acf_mol *c, char *f) {
  FILE *file = file_open(f,"r");
  mol_acf_fread(c,file);
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
