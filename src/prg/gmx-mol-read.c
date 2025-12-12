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
#include <cmn/message.h>
#include <cmn/pipe.h>
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/types.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Read molecule-type force-field data from gromacs topology file

   d    - molecular data structure
   line - storage for read line 
   f    - open file stream */
char* gmx_dat_read_moleculetype(struct gmx_mol *d, char *line, FILE *f) {
  char chr,name[80];
  /* read data */
  for (line=str_free(line),line=str_read_line_new(f); line;
       line=str_free(line),line=str_read_line_new(f)) {
    chr = str_char_first_nb(line,NULL);
    /* next keyword */
    if (chr=='[')
      break;
    /* non-commented line */
    if (chr && chr!=';') {
      if (sscanf(line,"%s%d",name,&(d->n_excl))!=2)
        msg_error("invalid format of molecular header",1);
      d->name = str_copy_new(name);
      }
    }
  return(line);
  }

/* Read atomic data from gromacs topology file

   d    - molecular data structure
   line - storage for read line 
   f    - open file stream */
char* gmx_dat_read_atoms(struct gmx_mol *d, char *line, FILE *f) {
  char chr,type[80],atom[80],res[80];
  unsigned i,j,r0,r1,ia=0;
  struct gmx_atom *a;
  struct gmx_res *r = NULL;
  struct queue *q,*p;
  if (!d)
    msg_error("'atoms' data block without preceeding 'moleculetype'",1);
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
      /* atom */
      a = gmx_atom_new(1);
      if (sscanf(line,"%*d%s%d%s%s%*d%lf%lf",type,
          &r1,res,atom,&(a->charge),&(a->mass))!=6) {
        if (sscanf(line,"%*d%s%d%s%s%*d%lf",type,
            &r1,res,atom,&(a->charge))!=5)
          msg_error("invalid format in 'atoms' data block",1);
        /* default mass value */
        a->mass = atom_mass_name_pdb(atom);
        }
      a->id = ia++;
      a->num = atom_num_pdb(atom);
      a->name = str_copy_new(atom);
      a->type = str_copy_new(type);
      /* residuum */
      if (!r || r0!=r1) {
        r = gmx_res_new(1);
        r->name = str_copy_new(res);
        p = queue_alloc();
        r->atom = (struct gmx_atom*)p;
        queue_add(q,r);
        }
      queue_add(p,a);
      r0 = r1;
      }
    }
  /* save data */
  d->n_resids = q->num;
  d->res = gmx_res_new(d->n_resids);
  for (i=0; i<d->n_resids; i++) {
    r = queue_get(q);
    p = (struct queue*)r->atom;
    r->n_atoms = p->num;
    r->atom = gmx_atom_new(r->n_atoms);
    for (j=0; j<r->n_atoms; j++) {
      a = queue_get(p);
      gmx_atom_copy(r->atom+j,a);
      gmx_atom_free(a,1);
      }
    queue_free(p);
    gmx_res_copy(d->res+i,r);
    gmx_res_free(r,1);
    }
  /* clean memory */
  queue_free(q);
  return(line);
  }

/* Read bonding data from gromacs topology file

   d    - storage for the bonding data
   nd   - number of bonding terms
   key  - data block keyword
   na   - number of atom IDs
   line - storage for read line 
   f    - open file stream */
char* gmx_dat_read_bonding(struct gmx_bond **d, unsigned *nd, char *key,
  unsigned na, char *line, FILE *f) {
  char **w,chr;
  unsigned i,n;
  struct gmx_bond *b,*t;
  struct queue *q;
  if (!d)
    msg_error_f("'%s' data block without preceeding 'moleculetype'",1,key);
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
        msg_error_f("inconsistent number of columns in '%s' data block",1,key);
      b = gmx_bond_new(1);
      /* atom IDs */
      b->n_atoms = na;
      b->ia = vec_ualloc(b->n_atoms);
      for (i=0; i<b->n_atoms; i++) {
        if (!type_read(w[i],&(b->ia[i]),TYPE_UINT))
          msg_error_f("cannot read atom ID in column #%d of '%s' data block",
            1,i+1,key);
        if (b->ia[i]<1)
          msg_error_f("invalid bonding ID in column #%d of '%s' data block",
            1,i+1,key);
        b->ia[i]--;
        }
      /* function type */
      if (!type_read(w[na],&(b->type),TYPE_SINT))
        msg_error_f("cannot read functional form ID in '%s' data block",1,key);
      /* explicit parameters */
      b->n_parms = n-na-1;
      b->parm = vec_falloc(b->n_parms);
      for (i=0; i<b->n_parms; i++)
        if (!type_read(w[na+i+1],&(b->parm[i]),TYPE_DOUBLE))
          msg_error_f("cannot read bonding parameter in '%s' data block",1,key);
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

/* Read list of exclusions from gromacs topology file

   d    - topology data structure
   line - storage for read line 
   f    - open file stream */
char* gmx_dat_read_exclusions(struct gmx_top *d, char *line, FILE *f) {
  char chr;
  unsigned i,nt;
  struct queue *q;
  struct gmx_bond *b,*t;
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
      /* list of exclusions */
      b = gmx_bond_new(1);
      b->ia = str_parse_uarray(line,&(b->n_atoms));
      for (i=0; i<b->n_atoms; i++) {
        if (b->ia[i]<1)
          msg_error_f("invalid bonding ID in column #%d of '%s' data block",
            1,i+1,"exclusions");
        b->ia[i]--;
        }
      queue_add(q,b);
      }
    }
  /* save data */
  if (q->num) {
    nt = d->n_exclusions + q->num;
    t = gmx_bond_new(nt);
    for (i=0; i<d->n_exclusions; i++)
      gmx_bond_copy(t+i,d->excl+i);
    for (i=d->n_exclusions; i<nt; i++) {
      b = queue_get(q);
      gmx_bond_copy(t+i,b);
      gmx_bond_free(b,1);
      }
    d->excl = gmx_bond_free(d->excl,d->n_exclusions);
    d->n_exclusions = nt;
    d->excl = t;
    }
  /* clean memory */
  queue_free(q);
  return(line);
  }

/* Read one molecule defintion from open gromacs topology file
 
   m    - molecular data structure
   line - storage for the read line
   f    - open file stream */
char* gmx_mol_fread_itp_one(struct gmx_mol *m, char *line, FILE *f) {
  char chr,*key;
  unsigned i;
  /* molecule type */
  if (!line)
    line = str_ffind_new(f,"[ moleculetype ]");
  line = gmx_dat_read_moleculetype(m,line,f);
  /* molecular data */
  while (line) {
    chr = str_char_first_nb(line,NULL);
    /* data block */
    if (chr=='[') {
      key = gmx_dat_read_top_key(line);
      /* atoms */
      if (str_compare(key,"atoms")) {
        line = gmx_dat_read_atoms(m,line,f);
        continue;
        }
      /* bonds */
      else if (str_compare(key,"bonds")) {
        line = gmx_dat_read_bonding(&(m->top->bond),&(m->top->n_bonds),
          "bonds",2,line,f);
        continue;
        }
      /* 1-4 pairs */
      else if (str_compare(key,"pairs")) {
        line = gmx_dat_read_bonding(&(m->top->pair),&(m->top->n_14_pairs),
          "pairs",2,line,f);
        continue;
        }
      /* angles */
      else if (str_compare(key,"angles")) {
        line = gmx_dat_read_bonding(&(m->top->angl),&(m->top->n_angles),
          "angles",3,line,f);
        continue;
        }
      /* dihedrals */
      else if (str_compare(key,"dihedrals")) {
        line = gmx_dat_read_bonding(&(m->top->dihe),&(m->top->n_dihedrals),
          "dihedrals",4,line,f);
        gmx_top_dihe_prop_impr(m->top);
        continue;
        }
      /* cmap */
      else if (str_compare(key,"cmap")) {
        line = gmx_dat_read_bonding(&(m->top->cmap),&(m->top->n_cmap),
          "cmap",5,line,f);
        continue;
        }
      /* constraints */
      else if (str_compare(key,"constraints")) {
        line = gmx_dat_read_bonding(&(m->top->cnst),&(m->top->n_constraints),
          "constraint",2,line,f);
        continue;
        }
      /* settles */
      else if (str_compare(key,"settles")) {
        line = gmx_dat_read_bonding(&(m->top->sttl),&(m->top->n_settles),
          "settle",1,line,f);
        continue;
        }
      /* exclusions */
      else if (str_compare(key,"exclusions")) {
        line = gmx_dat_read_exclusions(m->top,line,f);
        continue;
        }
      /* unknown keyword */
      else
        break;
      }
    line = str_free(line);
    line = str_read_line_new(f);
    }
  /* set number of atoms and residuum IDs */
  m->n_atoms = 0;
  for (i=0; i<m->n_resids; i++) {
    m->n_atoms += m->res[i].n_atoms;
    m->res[i].id = i;
    }
  /* set bonding IDs */
  gmx_mol_set_bonding(m);
  return(line);
  }

/* Read one molecule definition from gromacs topology file

   m      - molecular data structure
   name   - name of the file
   def    - array of pre-processor definitions
   n_defs - number of pre-processor definitions */
void gmx_mol_read_itp_one(struct gmx_mol *m, char *name,
  char **def, unsigned n_defs) {
  FILE *f;
  /* parse the file with preprocessor */
  f = gmx_read_pipe_open(name,def,n_defs);
  /* read data from file */
  gmx_mol_fread_itp_one(m,NULL,f);
  pipe_close(f);
  }

/* Read all molecule definitions from open gromacs topology file
 
   n - number of molecules (output)
   f - open file stream */
struct gmx_mol* gmx_mol_fread_itp(unsigned *n, FILE *f) {
  char chr,*line,*key;
  unsigned i;
  struct queue *q;
  struct gmx_mol *m,*t;
  /* read data from file */
  q = queue_alloc();
  line = str_read_line_new(f);
  while (line) {
    chr = str_char_first_nb(line,NULL);
    /* data block */
    if (chr=='[') {
      key = gmx_dat_read_top_key(line); 
      /* new molecule */
      if (str_compare(key,"moleculetype")) {
        m = gmx_mol_new(1);
        queue_add(q,m);
        line = gmx_mol_fread_itp_one(m,line,f);
        continue;
        }
      /* end of molecular data */
      else if (str_compare(key,"system")) {
        line = str_free(line);
        break;
        }
      /* unknown keyword */
      else
        msg_warn_f("unknown data block '%s' in gromacs topology",key);
      }
    line = str_free(line);
    line = str_read_line_new(f);
    }
  /* save data to array */
  (*n) = q->num;
  t = gmx_mol_new(*n);
  for (i=0; i<(*n); i++) {
    m = queue_get(q);
    gmx_mol_copy(t+i,m);
    gmx_mol_free(m,1);
    }
  /* clean memory */
  queue_free(q);
  return(t);
  }

/* Read all molecule definitions from gromacs topology file */
struct gmx_mol* gmx_mol_read_itp(char *name, unsigned *n, 
  char **def, unsigned n_defs) {
  FILE *f;
  struct gmx_mol *m;
  /* parse the file with preprocessor */
  f = gmx_read_pipe_open(name,def,n_defs);
  /* read data from file */
  m = gmx_mol_fread_itp(n,f);
  pipe_close(f);
  return(m);
  }

/* -------------------------------------------------------------------------- */
