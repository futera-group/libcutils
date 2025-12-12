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

#include <cmn/file.h>
#include <cmn/list.h>
#include <cmn/message.h>
#include <cmn/queue.h>
#include <cmn/string.h>
#include "mol/atom.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* assign value of internal coordinate from dictionary 
 
   t - dictionary with internal coordinates
   s - symbol of the coordinate
   v - value of the coordinate (output) */
void mol_zmt_assign_val(struct list *t, char *s, double *v) {
  double val = 0.0;
  short found = 0;
  struct ldata *p;
  struct zmt_icrd *d;
  for (p=t->first; p; p=p->l_next) {
    d = (struct zmt_icrd*)p->l_data;
    if (str_compare(d->sym,s)) {
      val = d->val;
      found = 1;
      break;
      }
    }
  if (!found)
    msg_error_f("value of internal coordinate \"%s\" not found",1,s);
  (*v) = val;
  }

/* read molecular data from zmt file
 
   z    - pointer to zmt molecular data struct 
   file - pointer to open of the file */
void mol_zmt_fread(struct zmt_mol *z, FILE *file) {
  char *line,sym[256],sym_b[256],sym_a[256],sym_d[256];
  unsigned id,nc = 0;
  struct queue *q;
  struct list *t;
  struct zmt_line *d;
  struct zmt_icrd *v;
  /* skip free lines */
  for (line=str_read_line_new(file); line;
       str_free(line),line=str_read_line_new(file)) {
    str_trim(line);
    if (line[0]!='\0' && line[0]!='\n')
      break;
    }
  /* read z-matrix */
  q = queue_alloc();
  z->n_atoms = 0;
  do {
    /* terminal blank line */
    str_trim(line);
    if (line[0]=='\0' || line[0]=='\n')
      break;
    /* read internal coordinates */
    d = mol_zmt_line_new();
    switch (z->n_atoms) {
      /* invalid format */
      case 0:
        if (sscanf(line,"%s",sym)!=1)
          msg_error("invalid format of line 1 in z-matrix",1);
        break;
      /* first atom */
      case 1:
        if (sscanf(line,"%s%u%s",sym,&(d->bond_id),sym_b)!=3)
          msg_error("invalid format of line 2 in z-matrix",1);
        d->bond_id--;
        d->bond_sym = str_copy_new(sym_b);
        nc = nc+1;
        break;
      /* second atom */
      case 2:
        if (sscanf(line,"%s%u%s%u%s",sym,&(d->bond_id),sym_b,
              &(d->angle_id),sym_a)!=5)
          msg_error("invalid format of line 3 in z-matrix",1);
        d->bond_id--;
        d->bond_sym = str_copy_new(sym_b);
        d->angle_id--;
        d->angle_sym = str_copy_new(sym_a);
        nc = nc+2;
        break;
      /* full specification */
      default:
        if (sscanf(line,"%s%u%s%u%s%u%s",sym,&(d->bond_id),sym_b,
              &(d->angle_id),sym_a,&(d->dihed_id),sym_d)!=7)
          msg_error_f("invalid format of line %d in z-matrix",1,z->n_atoms+1);
        d->bond_id--;
        d->bond_sym = str_copy_new(sym_b);
        d->angle_id--;
        d->angle_sym = str_copy_new(sym_a);
        d->dihed_id--;
        d->dihed_sym = str_copy_new(sym_d);
        nc = nc+3;
        break;
      }
    d->num = atom_num(sym);
    queue_add(q,d);
    z->n_atoms++;
    /* next line */
    str_free(line);
    line = str_read_line_new(file);
    }
  while (line);
  /* read internal coordinate values */
  t = list_alloc();
  for (line=str_read_line_new(file); line;
       str_free(line),line=str_read_line_new(file)) {
    /* blank line */
    str_trim(line);
    if (line[0]=='\0' || line[0]=='\n')
      break;
    /* coordinate value */
    v = mol_zmt_icrd_new();
    if (sscanf(line,"%s%lf",sym,&(v->val))!=2)
      msg_error("invalid format of internal coordinate values in z-matrix",1);
    v->sym = str_copy_new(sym);
    list_add_end_p(t,v);
    }
  /* convert queue to data struct */
  z->atom = mol_zmt_atom_new(z->n_atoms);
  id = 0;
  while (q->num) {
    d = queue_get(q);
    z->atom[id].num = d->num;
    z->atom[id].bond_id = d->bond_id;
    z->atom[id].angle_id = d->angle_id;
    z->atom[id].dihed_id = d->dihed_id;
    if (id>0)
      mol_zmt_assign_val(t,d->bond_sym,&(z->atom[id].bond_val));
    if (id>1)
      mol_zmt_assign_val(t,d->angle_sym,&(z->atom[id].angle_val));
    if (id>2)
      mol_zmt_assign_val(t,d->dihed_sym,&(z->atom[id].dihed_val));
    mol_zmt_line_free(d);
    id++;
    }
  /* clean memory */
  queue_free(q);
  list_free(t,NULL);
  }

/* read molecular data from zmt file
 
   z    - pointer to zmt molecular data struct 
   name - name of the file */
void mol_zmt_read(struct zmt_mol *z, char *name) {
  FILE *file;
  file = file_open(name,"r");
  mol_zmt_fread(z,file);
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
