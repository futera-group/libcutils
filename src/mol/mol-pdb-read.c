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

#include <string.h>
#include <cmn/file.h>
#include <cmn/list.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* first read of pdb file, save atom info in a queue for later reordering
 
   p    - pdb molecular data struct
   s    - list for data saving
   f    - open input filestream
   stop - specify if terminate program in case of error or not */
void mol_pdb_fread_unit(struct pdb_mol *p, struct list *s, 
 short stop, FILE *f) {
  unsigned length,resn,resn0,n,resid;
  char *line;
  struct pdb_unit *a,*last;
  /* initialization */
  resn = 0;
  resn0 = 0;
  resid = 1;
  n = 0;
  a = NULL;
  last = NULL;
  /* read file line by line */
  for (line=str_read_line_new(f); line;
       str_free(line),line=str_read_line_new(f)) {
    length = str_length(line);
    /* atomic data */
    if (str_sub_bfind(line,"ATOM") || str_sub_bfind(line,"HETATM")) {
      a = mol_pdb_unit_new();
      /* atom name */
      if (length<17) {
        if (stop)
          msg_error("incorrect pdb format - cannot read atom name",1);
        else
          return;
        }
      a->name = str_new(4);
      strncpy(a->name,line+12,4);
      a->name[4]='\0';
      /* residuum name */
      if (length<21) {
        if (stop)
          msg_error("incorrect pdb format - cannot read residuum name",1);
        else
          return;
        }
      a->resname = str_new(3);
      strncpy(a->resname,line+17,3);
      a->resname[3]='\0';
      /* residuum number */
      if (length<23 || sscanf(line+22,"%u",&resn)!=1) {
        if (stop)
          msg_error("incorrect pdb format - cannot read residuum number",1);
        else
          return;
        }
      /* coordinates */
      if (length<31 || sscanf(line+30,"%lf%lf%lf",
        &(a->coord[0]),&(a->coord[1]),&(a->coord[2]))!=3) {
        if (stop)
          msg_error("incorrect pdb format - cannot read coordinates",1);
        else
          return;
        }
      /* charge */
      if (length>54)
        sscanf(line+54,"%lf",&(a->charge));
      /* add atom to list */
      if (n>0 && (!str_compare(last->resname,a->resname) || resn!=resn0)) 
        resid++;
      resn0 = resn;
      a->resid = resid;
      list_add_end_p(s,a);
      last = a;
      n++;
      }
    /* termination flag */
    else if (str_sub_bfind(line,"TER")) 
      last->ter = 1;
    /* title */
    else if (str_sub_bfind(line,"TITLE") && length>6) {
      p->title = str_copy_new(line+6);
      str_trim(p->title);
      }
    /* end */
    else if (str_sub_bfind(line,"END"))
      break;
    }
  }

/* return number of different residui in the read units 
 
   s - list of the units */
unsigned mol_pdb_fread_nres(struct list *s) {
  unsigned i,n_max,n_res,*id;
  struct ldata *sp;
  /* initialization */
  n_max = 0;
  n_res = 0;
  /* largest residuum ID */
  sp = s->first;
  while (sp) {
    if (((struct pdb_unit*)(sp->l_data))->resid>n_max)
      n_max++;
    sp = sp->l_next;
    }
  /* count different residui */
  id = vec_ualloc(n_max+1);
  vec_uset(id,0,n_max+1);
  sp = s->first;
  while (sp) {
    id[((struct pdb_unit*)(sp->l_data))->resid]++;
    sp = sp->l_next;
    }
  for (i=0; i<(n_max+1); i++) {
    if (id[i])
      n_res++;
    }
  /* clean memory */
  vec_ufree(id);
  return(n_res);
  }

/* return number of atom in specific residuum
 
   s - list of the units
   r - residuum ID */
unsigned mol_pdb_fread_natoms(struct list *s, unsigned r) {
  unsigned n = 0;
  struct ldata *sp;
  sp = s->first;
  while (sp) {
    if (((struct pdb_unit*)(sp->l_data))->resid==r)
      n++;
    sp = sp->l_next;
    }
  return(n);
  }

/* read molecular data from open pdb file 
 
   p - pointer to pdb molecular data struct 
   t - specify if terminate program in case of error or not
   f - pointer to open PDB file */
int mol_pdb_fread(struct pdb_mol *p, short t, FILE *f) {
  unsigned ir,ia,it,ic;
  struct ldata *sp;
  struct list *s;
  struct pdb_unit *u;
  s = list_alloc();
  /* read data to list */
  mol_pdb_fread_unit(p,s,t,f);
  /* allocate residui */
  p->n_res = mol_pdb_fread_nres(s);
  p->res = mol_pdb_res_new(p->n_res);
  /* allocate atoms */
  it = 1;
  for (ir=0; ir<p->n_res; ir++) {
    p->res[ir].n_atoms = mol_pdb_fread_natoms(s,ir+1);
    p->res[ir].atom = mol_pdb_atom_new(p->res[ir].n_atoms);
    /* copy data */
    ia = 0;
    for (sp=s->first; sp; sp = sp->l_next) {
      u = (struct pdb_unit*)(sp->l_data);
      if (u->resid==(ir+1)) {
        /* residuum data */
        p->res[ir].name = u->resname;
        u->resname = NULL;
        p->res[ir].id = u->resid;
        if (u->ter)
          p->res[ir].ter = u->ter;
        /* atom data */
        for (ic=0; ic<3; ic++)
          p->res[ir].atom[ia].coord[ic] = u->coord[ic];
        p->res[ir].atom[ia].name = u->name;
        u->name = NULL;
        p->res[ir].atom[ia].charge = u->charge;
        p->res[ir].atom[ia].id = it;
        ia++;
        it++;
        }
      }
    }
  /* clean memory */
  list_free(s,NULL);
  return(1);
  }

/* read molecular data from pdb file 
 
   p    - pointer to pdb molecular data struct 
   name - name of the file */
void mol_pdb_read(struct pdb_mol *p, char *name) {
  FILE *file;
  file = file_open(name,"r");
  mol_pdb_fread(p,1,file);
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
