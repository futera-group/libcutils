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
#include <cmn/pipe.h>
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/types.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Read data-block keyword from line extracted from Gromacs file 
 
   line - line from the file */
char *gmx_dat_read_top_key(char *line) {
  static char key[80];
  str_copy(key,line);
  str_ltrim_mark(key,'[');
  str_rtrim_mark(key,']');
  str_trim(key);
  return(key);
  }

/* -------------------------------------------------------------------------- */

/* Read default force-field setting from gromacs topology file

   d    - force-field data structure
   line - storage for read line 
   f    - open file stream */
char* gmx_dat_read_defaults(struct gmx_ff *d, char *line, FILE *f) {
  char chr,gen[80];
  /* read data */
  for (line=str_free(line),line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    chr = str_char_first_nb(line,NULL);
    /* next keyword */
    if (chr=='[')
      break;
    /* non-commented line */
    if (chr && chr!=';') {
      if (sscanf(line,"%hd%hd%s%lf%lf",&(d->nb_func),&(d->comb_rule),
          gen,&(d->scale_lj),&(d->scale_qq))!=5)
        msg_error("invalid format in 'defaults' data block",1);
      if (str_compare(gen,"yes"))
        d->gen_pairs = 1;
      else if (str_compare(gen,"no"))
        d->gen_pairs = 0;
      else
        msg_error("invalid 'gen-pairs' option value \"%s\"",1);
      }
    }
  return(line);
  }

/* Read atom-type force-field data from gromacs topology file

   d    - force-field data structure
   line - storage for read line 
   f    - open file stream */
char* gmx_dat_read_atomtypes(struct gmx_ff *d, char *line, FILE *f) {
  char chr,sym[80],type[80];
  unsigned i,n;
  struct gmx_atom *a,*t;
  struct queue *q;
  q = queue_alloc();
  /* read data */
  for (line=str_free(line),line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    chr = str_char_first_nb(line,NULL);
    /* next keyword */
    if (chr=='[')
      break;
    /* non-commented line */
    if (chr && chr!=';') {
      a = gmx_atom_new(1);
      if (sscanf(line,"%s%u%lf%lf%s%lf%lf",type,&(a->num),&(a->mass),
          &(a->charge),sym,&(a->sigma),&(a->epsilon))!=7)
        msg_error("invalid format in 'atomtypes' data block",1);
      a->particle = gmx_ff_particle_type_id(sym);
      a->type = str_copy_new(type);
      queue_add(q,a);
      }
    }
  /* save data */
  if (q->num) {
    n = d->dat->n_atom_types + q->num;
    t = gmx_atom_new(n);
    for (i=0; i<d->dat->n_atom_types; i++)
      gmx_atom_copy(t+i,d->dat->type+i);
    for (i=d->dat->n_atom_types; i<n; i++) {
      a = queue_get(q);
      gmx_atom_copy(t+i,a);
      gmx_atom_free(a,1);
      }
    d->dat->type = gmx_atom_free(d->dat->type,d->dat->n_atom_types);
    d->dat->n_atom_types = n;
    d->dat->type = t;
    }
  /* clean memory */
  queue_free(q);
  return(line);
  }

/* Read bonding-term force-field data from gromacs topology file

   d    - storage for the data
   nd   - length of the data array
   key  - data file keyword
   na   - number of atom types
   line - storage for read line 
   f    - open file stream */
char* gmx_dat_read_bondtypes(struct gmx_bond **d, unsigned *nd, 
  char *key, unsigned na, char *line, FILE *f) {
  char **w,chr;
  unsigned i,n;
  struct gmx_bond *b,*t;
  struct queue *q;
  q = queue_alloc();
  /* read data */
  for (line=str_free(line),line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    chr = str_char_first_nb(line,NULL);
    /* next keyword */
    if (chr=='[')
      break;
    /* non-commented line */
    if (chr && chr!=';') {
      str_rtrim_mark(line,';');
      w = str_split_s(line,"\t\n ",&n);
      if (n<(na+1))
        msg_error_f("inconsistent number of columns in '%s'"
          " data block",1,key);
      b = gmx_bond_new(1);
      b->n_atoms = na;
      b->name = vec_talloc(sizeof(char*),b->n_atoms);
      for (i=0; i<na; i++)
        b->name[i] = str_copy_new(w[i]);
      if (!type_read(w[na],&(b->type),TYPE_SINT))
        msg_error_f("cannot read bonding functional type in"
          " '%s' data block",1,key);
      b->n_parms = n - b->n_atoms - 1;
      b->parm = vec_falloc(b->n_parms);
      for (i=0; i<b->n_parms; i++)
        if (!type_read(w[na+i+1],&(b->parm[i]),TYPE_DOUBLE))
          msg_error_f("cannot read bonding parameters #%d in"
            " '%s' data block",1,i+1,key);
      queue_add(q,b);
      vec_sfree(w,n);
      }
    }
  /* save data */
  if (q->num) {
    n = (*nd) + q->num;
    t = gmx_bond_new(n);
    for (i=0; i<(*nd); i++)
      gmx_bond_copy(t+i,(*d)+i);
    for (i=(*nd); i<n; i++) {
      b = queue_get(q);
      gmx_bond_copy(t+i,b);
      gmx_bond_free(b,1);
      }
    (*d) = gmx_bond_free(*d,*nd);
    (*nd) = n;
    (*d) = t;
    }
  /* clean memory */
  queue_free(q);
  return(line);
  }

/* Read implicit solvent parameters from gromacs topology file

   d    - storage for the data
   line - storage for read line 
   f    - open file stream */
char* gmx_dat_read_igb(struct gmx_ff *d, char *line, FILE *f) {
  char chr,type[80];
  double p1,p2,p3,p4,p5;
  unsigned i,n;
  struct gmx_bond *b,*t;
  struct queue *q;
  q = queue_alloc();
  /* read data */
  for (line=str_free(line),line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    chr = str_char_first_nb(line,NULL);
    /* next keyword */
    if (chr=='[')
      break;
    /* non-commented line */
    if (chr && chr!=';') {
      if (sscanf(line,"%s%lf%lf%lf%lf%lf",type,&p1,&p2,&p3,&p4,&p5)!=6)
        msg_error("invalid format in 'implicit_genborn_params' data block",1);
      b = gmx_bond_new(1);
      b->n_atoms = 1;
      b->name = vec_talloc(sizeof(char*),1);
      b->name[0] = str_copy_new(type);
      b->n_parms = 5;
      b->parm = vec_falloc(b->n_parms);
      b->parm[0] = p1;
      b->parm[1] = p2;
      b->parm[2] = p3;
      b->parm[3] = p4;
      b->parm[4] = p5;
      queue_add(q,b);
      }
    }
  /* save data */
  if (q->num) {
    n = d->dat->n_igb_parms + q->num;
    t = gmx_bond_new(n);
    for (i=0; i<d->dat->n_igb_parms; i++)
      gmx_bond_copy(t+i,d->dat->igbp+i);
    for (i=d->dat->n_igb_parms; i<n; i++) {
      b = queue_get(q);
      gmx_bond_copy(t+i,b);
      gmx_bond_free(b,1);
      }
    d->dat->igbp = gmx_bond_free(d->dat->igbp,d->dat->n_igb_parms);
    d->dat->n_igb_parms = n;
    d->dat->igbp = t;
    }
  /* clean memory */
  queue_free(q);
  return(line);
  }

/* Read system data block from gromacs topology file

   d    - gromacs data structure
   line - storage for read line 
   f    - open file stream */
char* gmx_dat_read_system(struct gmx_dat *d, char *line, FILE *f) {
  char chr;
  /* read data */
  for (line=str_free(line),line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    chr = str_char_first_nb(line,NULL);
    /* next keyword */
    if (chr=='[')
      break;
    /* non-commented line */
    if (chr && chr!=';') {
      str_rtrim_mark(line,';');
      str_trim(line);
      d->title = str_copy_new(line);
      }
    }
  return(line);
  }

/* Read molecular array definition from gromacs topology file

   q    - storage for the molecular records
   line - storage for read line 
   f    - open file stream */
char* gmx_dat_read_molecules(struct queue *q, char *line, FILE *f) {
  char chr,name[80];
  struct gmx_frag *m;
  /* read data */
  for (line=str_free(line),line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    chr = str_char_first_nb(line,NULL);
    /* next keyword */
    if (chr=='[')
      break;
    /* non-commented line */
    if (chr && chr!=';') {
      m = gmx_frag_new(1);
      if (sscanf(line,"%s%u",name,&(m->n_rep))!=2)
        msg_error("invalid format in 'molecules' data block",1);
      m->name = str_copy_new(name);
      queue_add(q,m);
      }
    }
  return(line);
  }

/* Read data blocks from open gromacs topology file

   d - gromacs data structure
   q - storage for molecular topology
   t - storage for molecular array
   f - open file stream */
void gmx_dat_read_top_data(struct gmx_dat *d, struct queue *q, 
  struct queue *t, FILE *f) {
  char chr,*line,*key;
  struct gmx_mol *m = NULL;
  line = str_read_line_new(f);
  while (line) {
    chr = str_char_first_nb(line,NULL);
    /* data block */
    if (chr=='[') {
      key = gmx_dat_read_top_key(line);
      /* default setting */
      if (str_compare(key,"defaults")) {
        line = gmx_dat_read_defaults(d->ff,line,f);
        continue;
        }
      /* atom types */
      else if (str_compare(key,"atomtypes")) {
        line = gmx_dat_read_atomtypes(d->ff,line,f);
        continue;
        }
      /* 1-4 pair types */
      else if (str_compare(key,"pairtypes")) {
        line = gmx_dat_read_bondtypes(&(d->ff->dat->pair),
          &(d->ff->dat->n_14_pairs),key,2,line,f);
        continue;
        }
      /* non-bonding parameters */
      else if (str_compare(key,"nonbond_params")) {
        line = gmx_dat_read_bondtypes(&(d->ff->dat->nbpr),
          &(d->ff->dat->n_nb_pairs),key,2,line,f);
        continue;
        }
      /* bond types */
      else if (str_compare(key,"bondtypes")) {
        line = gmx_dat_read_bondtypes(&(d->ff->dat->bond),
          &(d->ff->dat->n_bonds),key,2,line,f);
        continue;
        }
      /* constrain types */
      else if (str_compare(key,"constrainttypes")) {
        line = gmx_dat_read_bondtypes(&(d->ff->dat->cnst),
          &(d->ff->dat->n_constraints),key,2,line,f);
        continue;
        }
      /* angle types */
      else if (str_compare(key,"angletypes")) {
        line = gmx_dat_read_bondtypes(&(d->ff->dat->angl),
          &(d->ff->dat->n_angles),key,3,line,f);
        continue;
        }
      /* dihedral types */
      else if (str_compare(key,"dihedraltypes")) {
        line = gmx_dat_read_bondtypes(&(d->ff->dat->dihe),
          &(d->ff->dat->n_dihedrals),key,4,line,f);
        gmx_top_dihe_prop_impr(d->ff->dat);
        continue;
        }
      /* implicit solvent parameters */
      else if (str_compare(key,"implicit_genborn_params")) {
        line = gmx_dat_read_igb(d->ff,line,f);
        continue;
        }
      /* cmap types */
      else if (str_compare(key,"cmaptypes")) {
        line = gmx_dat_read_bondtypes(&(d->ff->dat->cmap),
          &(d->ff->dat->n_cmap),key,5,line,f);
        continue;
        }
      /* molecule types */
      else if (str_compare(key,"moleculetype")) {
        m = gmx_mol_new(1);
        queue_add(q,m);
        line = gmx_mol_fread_itp_one(m,line,f);
        continue;
        }
      /* system */
      else if (str_compare(key,"system")) {
        line = gmx_dat_read_system(d,line,f);
        continue;
        }
      /* molecules */
      else if (str_compare(key,"molecules")) {
        line = gmx_dat_read_molecules(t,line,f);
        continue;
        }
      /* unknown keyword */
      else
        msg_warn_f("unknown data block '%s' in gromacs topology",key);
      }
    line = str_free(line);
    line = str_read_line_new(f);
    }
  }

/* Read data from gromacs topology file

   d    - gromacs data structure
   name - name of the file */
void gmx_dat_read_top(struct gmx_dat *d, char *name) {
  unsigned i;
  FILE *f;
  struct gmx_mol *m;
  struct gmx_frag *s;
  struct queue *q,*t;
  /* parse the file with preprocessor */
  f = gmx_read_pipe_open(name,d->ff->def,d->ff->n_defs);
  /* read data from file */
  q = queue_alloc();
  t = queue_alloc();
  gmx_dat_read_top_data(d,q,t,f);
  pipe_close(f);
  /* save data */
  d->n_mols = q->num;
  d->mol = gmx_mol_new(d->n_mols);
  for (i=0; i<d->n_mols; i++) {
    m = queue_get(q);
    gmx_mol_copy(d->mol+i,m);
    gmx_mol_free(m,1);
    }
  /* system definition */
  d->n_frags = t->num;
  d->frag = gmx_frag_new(d->n_frags);
  for (i=0; i<d->n_frags; i++) {
    s = queue_get(t);
    gmx_frag_copy(d->frag+i,s);
    d->frag[i].mol = gmx_dat_mol_get(d,d->frag[i].name);
    gmx_frag_free(s,1);
    }
  gmx_dat_update_nums(d);
  /* clean memory */
  queue_free(q);
  queue_free(t);
  }

/* -------------------------------------------------------------------------- */
