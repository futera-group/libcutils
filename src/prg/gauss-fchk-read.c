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
#include <string.h>
#include <cmn/file.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/types.h>
#include <cmn/vector.h>
#include <qmc/basis.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* find and read one single item in gaussian formatted checkpoint file
 
   key  - keyword that is expected in file
   t    - type of read data
   d    - pointer to data storage
   name - name of the file
   file - pointer to open file */
void gauss_fchk_read_item(char *key, short t, void *d, char *name, FILE *file) {
  char *line;
  short pass = 0;
  line = str_ffind_b_new(file,key);
  while (!line) {
    if (pass)
      msg_error_f("cannot find \"%s\" in %s file",1,key,name);
    rewind(file);
    line = str_free(line);
    line = str_ffind_b_new(file,key);
    pass = 1;
    }
  switch (t) {
    case TYPE_INT:
     if (sscanf(line+50,"%d",(int*)d)!=1)
       msg_error_f("cannot read \"%s\" value in %s file",1,key,name);
     break;
    case TYPE_UINT:
     if (sscanf(line+50,"%u",(unsigned*)d)!=1)
       msg_error_f("cannot read \"%s\" value in %s file",1,key,name);
     break;
    default:
      msg_error_f("unsupported type in gauss_fchk_read_item"
        " (key = %s, t = %d, name = %s)",1,key,t,name);
    }
  line = str_free(line);
  }

/* find and read vector item in gaussian formatted checkpoint file
 
   key  - keyword that is expected in file
   t    - type of read data
   n    - lenght of data array (output)
   ne   - expectation of data array lenght (ignored if 0)
   wt   - specify if the data are mandatory (wanted = 1)
   name - name of the file
   file - pointer to open file */
void* gauss_fchk_read_item_v(char *key, short t, unsigned *n, unsigned ne,
  short wt, char *name, FILE *file) {
  char *line;
  unsigned i,j,nb,nn;
  short pass = 0;
  void *d = NULL;
  line = str_ffind_b_new(file,key);
  while (!line) {
    if (pass) {
      if (wt)
        msg_error_f("cannot find \"%s\" in %s file",1,key,name);
      else
        return(NULL);
      }
    rewind(file);
    line = str_free(line);
    line = str_ffind_b_new(file,key);
    pass = 1;
    }
  if (sscanf(line+50,"%u",n)!=1)
    msg_error_f("cannot read \"%s\" value in %s file",1,key,name);
  switch (t) {
    case TYPE_UINT:
      d = vec_ualloc(*n);
      nb = ((*n)%6 ? (*n)/6+1 : (*n)/6);
      line = str_free(line);
      for (i=0; i<nb; i++) {
        line = str_read_line_new(file);
        if (!line)
          msg_error_f("unexpected end of %s file while \"%s\"",1,name,key);
        nn = (i==(nb-1) ? ((*n)%6 ? (*n)%6 : 6) : 6);
        for (j=0; j<nn; j++)
          if (sscanf(line+12*j,"%u",&(((unsigned*)d)[6*i+j]))!=1)
            msg_error_f("cannot read \"%s\" data in %s file",1,key,name);
        line = str_free(line);
        }
      break;
    case TYPE_SINT:
      d = vec_sialloc(*n);
      nb = ((*n)%6 ? (*n)/6+1 : (*n)/6);
      line = str_free(line);
      for (i=0; i<nb; i++) {
        line = str_read_line_new(file);
        if (!line)
          msg_error_f("unexpected end of %s file while \"%s\"",1,name,key);
        nn = (i==(nb-1) ? ((*n)%6 ? (*n)%6 : 6) : 6);
        for (j=0; j<nn; j++)
          if (sscanf(line+12*j,"%hd",&(((short*)d)[6*i+j]))!=1)
            msg_error_f("cannot read \"%s\" data in %s file",1,key,name);
        line = str_free(line);
        }
      break;
    case TYPE_DOUBLE:
      d = vec_falloc(*n);
      nb = ((*n)%5 ? (*n)/5+1 : (*n)/5);
      line = str_free(line);
      for (i=0; i<nb; i++) {
        line = str_read_line_new(file);
        if (!line)
          msg_error_f("unexpected end of %s file while \"%s\"",1,name,key);
        nn = (i==(nb-1) ? ((*n)%5 ? (*n)%5 : 5) : 5);
        for (j=0; j<nn; j++)
          if (sscanf(line+16*j,"%lf",&(((double*)d)[5*i+j]))!=1)
            msg_error_f("cannot read \"%s\" data in %s file",1,key,name);
        line = str_free(line);
        }
      break;
    case TYPE_STRING:
      d = str_new(*n);
      nb = ((*n)%5 ? (*n)/5+1 : (*n)/5);
      line = str_free(line);
      for (i=0; i<nb; i++) {
        line = str_read_line_new(file);
        if (!line)
          msg_error_f("unexpected end of %s file while \"%s\"",1,name,key);
        nn = (i==(nb-1) ? ((*n)%5 ? (*n)%5 : 5) : 5);
        for (j=0; j<nn; j++)
          strncpy(((char**)d)[5*i+j],line+12*j,13);
        line = str_free(line);
        }
      break;
    default:
      msg_error_f("unsupported type in gauss_fchk_read_item_v (id = %d)",1,t);
    }
  if (ne && ((*n)!=ne))
    msg_error_f("inconsistent lenght of \"%s\" array"
     " (%d found / %d expected)",1,key,*n,ne);
  return(d);
  }

/* find out if gaussian data corresponds to restricted electronic system
 
   d - pointer to gaussian data struct */
short gauss_fchk_mthd_is_r(struct gauss_dat *d) {
  short rr = 0;
  if (d->job_mthd[0]=='R')
    rr = 1;
  else if (d->job_mthd[0]=='U')
    rr = 0;
  else
    msg_error_f("cannot decide type of method from \"%s\" flag",1,d->job_mthd);
  return(rr);
  }

/* read basis set data from gaussian formatted checkpoint file
 
   d    - pointer to gaussian data struct
   name - name of the file
   file - pointer to open file stream */
void gauss_fchk_read_basis(struct gauss_dat *d, char *name, FILE *file) {
  double *prim_exp,*cont_coeff,*cont_coeff_sp,*shell_coord;
  unsigned i,j,k,n,id_s,id_x,*prim_num,*atom_map;
  short *shell_type;
  struct basis *b;
  struct basis_center *c;
  struct basis_shell *s;
  rewind(file);
  /* clean basis set struct */
  d->bs = basis_free(d->bs);
  d->bs = basis_new();
  b = d->bs;
  /* read basis set specifiers */
  gauss_fchk_read_item("Number of basis functions",TYPE_UINT,
    &(b->n_bfce),name,file);
  gauss_fchk_read_item("Number of independent functions",TYPE_UINT,
    &(b->n_ibfce),name,file);
  gauss_fchk_read_item("Number of contracted shells",TYPE_UINT,
    &(b->n_cont_shells),name,file);
  gauss_fchk_read_item("Number of primitive shells",TYPE_UINT,
    &(b->n_prim_shells),name,file);
  gauss_fchk_read_item("Pure/Cartesian d shells",TYPE_UINT,
    &(b->n_pure_d),name,file);
  gauss_fchk_read_item("Pure/Cartesian f shells",TYPE_UINT,
    &(b->n_pure_f),name,file);
  gauss_fchk_read_item("Highest angular momentum",TYPE_UINT,
    &(b->max_ang_mom),name,file);
  gauss_fchk_read_item("Largest degree of contraction",TYPE_UINT,
    &(b->max_bfce_cont),name,file);
  /* read basis function data */
  shell_type = gauss_fchk_read_item_v("Shell types",TYPE_SINT,
    &n,b->n_cont_shells,1,name,file);
  prim_num = gauss_fchk_read_item_v("Number of primitives per shell",TYPE_UINT,
    &n,b->n_cont_shells,1,name,file);
  atom_map = gauss_fchk_read_item_v("Shell to atom map",TYPE_UINT,
    &n,b->n_cont_shells,1,name,file);
  prim_exp = gauss_fchk_read_item_v("Primitive exponents",TYPE_DOUBLE,
    &n,b->n_prim_shells,1,name,file);
  cont_coeff = gauss_fchk_read_item_v("Contraction coefficients",TYPE_DOUBLE,
    &n,b->n_prim_shells,1,name,file);
  cont_coeff_sp = gauss_fchk_read_item_v("P(S=P) Contraction coefficients",
    TYPE_DOUBLE,&n,b->n_prim_shells,0,name,file);
  shell_coord = gauss_fchk_read_item_v("Coordinates of each shell",TYPE_DOUBLE,
    &n,3*b->n_cont_shells,1,name,file);
  /* convert arrays to basis data structs */
  id_s = id_x = 0;
  b->n_centers = d->n_atoms;
  b->center = basis_center_new(b->n_centers);
  for (i=0; i<b->n_centers; i++) {
    c = b->center+i;
    /* atomic number and coordinates */
    c->type = d->atom[i].num;
    for (j=0; j<3; j++)
      c->coord[j] = d->atom[i].coord[j];
    /* number of shells */
    j = id_s;
    c->n_shells = 0;
    while (j<b->n_cont_shells && atom_map[j]==atom_map[id_s]) {
      c->n_shells++;
      j++;
      }
    /* atomic shell data */
    c->shell = basis_shell_new(c->n_shells);
    for (j=0; j<c->n_shells; j++) {
      s = c->shell+j;
      s->type = gauss_bs_shell_type(shell_type[id_s]);
      s->n_prim = prim_num[id_s];
      s->exp = vec_falloc(s->n_prim);
      s->cf1 = vec_falloc(s->n_prim);
      if (s->type==BASIS_SHELL_SP)
        s->cf2 = vec_falloc(s->n_prim);
      for (k=0; k<s->n_prim; k++) {
        s->exp[k] = prim_exp[id_x];
        s->cf1[k] = cont_coeff[id_x];
        if (s->type==BASIS_SHELL_SP)
          s->cf2[k] = cont_coeff_sp[id_x];
        id_x++;
        }
      id_s++;
      }
    }
  basis_def_parm(b);
  /* clean memory */
  vec_sifree(shell_type);
  vec_ufree(prim_num);
  vec_ufree(atom_map);
  vec_ffree(prim_exp);
  vec_ffree(cont_coeff);
  vec_ffree(cont_coeff_sp);
  vec_ffree(shell_coord);
  }

/* read MOs from gaussian formatted checkpoint file
 
   d    - pointer to gaussian data struct
   name - name of the file
   file - pointer to open file stream */
void gauss_fchk_read_mo(struct gauss_dat *d, char *name, FILE *file) {
  double *energy,*coeff;
  unsigned i,n;
  struct basis *b = d->bs;
  /* alpha MOs */
  energy = gauss_fchk_read_item_v("Alpha Orbital Energies",TYPE_DOUBLE,
    &n,b->n_ibfce,1,name,file);
  coeff = gauss_fchk_read_item_v("Alpha MO coefficients",TYPE_DOUBLE,
    &n,b->n_ibfce*b->n_bfce,1,name,file);
  d->mo_a = vec_talloc(b->n_ibfce,sizeof(struct gauss_mo));
  for (i=0; i<b->n_ibfce; i++) {
    d->mo_a[i].energy = energy[i];
    d->mo_a[i].coeff = vec_fcopy_new(coeff+i*b->n_bfce,b->n_bfce);
    }
  energy = vec_ffree(energy);
  coeff = vec_ffree(coeff);
  /* beta MOs */
  if (!gauss_fchk_mthd_is_r(d)) {
    energy = gauss_fchk_read_item_v("Beta Orbital Energies",TYPE_DOUBLE,
      &n,b->n_ibfce,1,name,file);
    coeff = gauss_fchk_read_item_v("Beta MO coefficients",TYPE_DOUBLE,
      &n,b->n_ibfce*b->n_bfce,1,name,file);
    d->mo_b = vec_talloc(b->n_ibfce,sizeof(struct gauss_mo));
    for (i=0; i<b->n_ibfce; i++) {
      d->mo_b[i].energy = energy[i];
      d->mo_b[i].coeff = vec_fcopy_new(coeff+i*b->n_bfce,b->n_bfce);
      }
    energy = vec_ffree(energy);
    coeff = vec_ffree(coeff);
    }
  }

/* read atomic data from gaussian formatted checkpoint file
 
   d    - pointer to gaussian data struct
   name - name of the file
   file - pointer to open file stream */
void gauss_fchk_read_atom(struct gauss_dat *d, char *name, FILE *file) {
  double *coord=NULL,*coord_fix=NULL,*grad=NULL,*ch_mul=NULL,*ch_nuc=NULL;
  double *ch_nuc_e=NULL,*ch_mm=NULL,*ch_mm_m=NULL,*mass=NULL,*nuc_qmom=NULL,*nuc_gfac=NULL;
  unsigned *atom_num=NULL,*nuc_spin=NULL,*res_num=NULL,*res_info=NULL;
  unsigned *frag_info=NULL,*weight=NULL,*type_id=NULL,*type_id_m=NULL;
  char **type_name=NULL,**type_name_m=NULL;
  unsigned i,j,n;
  /* allocate memory */
  d->atom = vec_talloc(d->n_atoms,sizeof(struct gauss_atom));
  /* read atomic data */
  atom_num = gauss_fchk_read_item_v("Atomic numbers",TYPE_UINT,
    &n,d->n_atoms,1,name,file);
  coord = gauss_fchk_read_item_v("Current cartesian coordinates",TYPE_DOUBLE,
    &n,3*d->n_atoms,1,name,file);
  coord_fix = gauss_fchk_read_item_v("Constraint Structure",TYPE_DOUBLE,
    &n,3*d->n_atoms,0,name,file);
  grad = gauss_fchk_read_item_v("Cartesian Gradient",TYPE_DOUBLE,
    &n,3*d->n_atoms,1,name,file);
  ch_nuc = gauss_fchk_read_item_v("Nuclear charges",TYPE_DOUBLE,
    &n,d->n_atoms,1,name,file);
  ch_nuc_e = gauss_fchk_read_item_v("Nuclear ZEff",TYPE_DOUBLE,
    &n,d->n_atoms,1,name,file);
  ch_mm = gauss_fchk_read_item_v("MM charges",TYPE_DOUBLE,
    &n,d->n_atoms,1,name,file);
  ch_mm_m = gauss_fchk_read_item_v("Atom Modified MM Charges",TYPE_DOUBLE,
    &n,d->n_atoms,1,name,file);
  ch_mul = gauss_fchk_read_item_v("Mulliken Charges",TYPE_DOUBLE,
    &n,d->n_atoms,0,name,file);
  mass = gauss_fchk_read_item_v("Real atomic weights",TYPE_DOUBLE,
    &n,d->n_atoms,1,name,file);
  type_id = gauss_fchk_read_item_v("Int Atom Types",TYPE_UINT,
    &n,d->n_atoms,1,name,file);
  type_id_m = gauss_fchk_read_item_v("Int Atom Modified Types",TYPE_UINT,
    &n,d->n_atoms,1,name,file);
  weight = gauss_fchk_read_item_v("Integer atomic weights",TYPE_UINT,
    &n,d->n_atoms,1,name,file);
  res_info = gauss_fchk_read_item_v("Atom residue info",TYPE_UINT,
    &n,d->n_atoms,0,name,file);
  frag_info = gauss_fchk_read_item_v("Atom fragment info",TYPE_UINT,
    &n,d->n_atoms,1,name,file);
  res_num = gauss_fchk_read_item_v("Atom residue num",TYPE_UINT,
    &n,d->n_atoms,1,name,file);
  nuc_spin = gauss_fchk_read_item_v("Nuclear spins",TYPE_UINT,
    &n,d->n_atoms,1,name,file);
  nuc_qmom = gauss_fchk_read_item_v("Nuclear QMom",TYPE_DOUBLE,
    &n,d->n_atoms,1,name,file);
  nuc_gfac = gauss_fchk_read_item_v("Nuclear GFac",TYPE_DOUBLE,
    &n,d->n_atoms,1,name,file);
  for (i=0; i<d->n_atoms; i++) {
    if (atom_num)
      d->atom[i].num = atom_num[i];
    if (ch_nuc)
      d->atom[i].charge[GAUSS_CH_NC] = ch_nuc[i];
    if (ch_nuc_e)
      d->atom[i].charge[GAUSS_CH_NEFF] = ch_nuc_e[i];
    if (ch_mm)
      d->atom[i].charge[GAUSS_CH_MM] = ch_mm[i];
    if (ch_mm_m)
      d->atom[i].charge[GAUSS_CH_MM_MOD] = ch_mm_m[i];
    if (ch_mul)
      d->atom[i].charge[GAUSS_CH_MULL] = ch_mul[i];
    if (mass)
      d->atom[i].mass = mass[i];
    if (type_id)
      d->atom[i].type_id[0] = type_id[i];
    if (type_id_m)
      d->atom[i].type_id[1] = type_id_m[i];
    if (weight)
      d->atom[i].weight = weight[i];
    d->atom[i].res_info = (res_info ? res_info[i] : 0);
    if (res_num)
      d->atom[i].res_num = res_num[i];
    if (frag_info)
      d->atom[i].frag_info = frag_info[i];
    if (nuc_spin)
      d->atom[i].nuc_spin = nuc_spin[i];
    if (nuc_qmom)
      d->atom[i].nuc_qmom = nuc_qmom[i];
    if (nuc_gfac)
      d->atom[i].nuc_gfac = nuc_gfac[i];
    if (type_name)
      sprintf(d->atom[i].type_name[0],"%s",type_name[i]);
    if (type_name_m)
      sprintf(d->atom[i].type_name[1],"%s",type_name_m[i]);
    if (coord)
      for (j=0; j<3; j++)
        d->atom[i].coord[j] = coord[3*i+j];
    if (coord_fix)
      for (j=0; j<3; j++)
        d->atom[i].coord_fix[j] = coord_fix[3*i+j];
    if (grad)
      for (j=0; j<3; j++)
        d->atom[i].grad[j] = grad[3*i+j];
    }
  /* clean memory */
  coord = vec_ffree(coord);
  coord_fix = vec_ffree(coord_fix);
  grad = vec_ffree(grad);
  atom_num = vec_ufree(atom_num);
  ch_nuc = vec_ffree(ch_nuc);
  ch_nuc_e = vec_ffree(ch_nuc_e);
  ch_mm = vec_ffree(ch_mm);
  ch_mm_m = vec_ffree(ch_mm_m);
  ch_mul = vec_ffree(ch_mul);
  mass = vec_ffree(mass);
  type_id = vec_ufree(type_id);
  type_id_m = vec_ufree(type_id_m);
  weight = vec_ufree(weight);
  res_info = vec_ufree(res_info);
  res_num = vec_ufree(res_num);
  frag_info = vec_ufree(frag_info);
  nuc_spin = vec_ufree(nuc_spin);
  nuc_qmom = vec_ffree(nuc_qmom);
  nuc_gfac = vec_ffree(nuc_gfac);
  type_name = vec_sfree(type_name,d->n_atoms);
  type_name_m = vec_sfree(type_name_m,d->n_atoms);
  }

/* read data from gaussian formatted checkpoint file
 
   d - pointer to gaussian data struct
   f - file name */
void gauss_fchk_read(struct gauss_dat *d, char *f) {
  char *line,sym[80],mthd[80],basis[80];
  FILE *file;
  file = file_open(f,"r");
  /* header */
  d->job_title = str_read_line_new(file);
  if (!d->job_title)
    msg_error_f("gaussian formatted checkpoint file %s is empty",1,f);
  str_trim(d->job_title);
  line = str_read_line_new(file);
  if (!line)
    msg_error_f("unexpected end of %s file while reading job type",1,f);
  if (sscanf(line,"%s%s%s",sym,mthd,basis)!=3)
    msg_error_f("invalid format of job type specification in %s",1,f);
  d->job_mthd = str_copy_new(mthd);
  d->job_basis = str_copy_new(basis);
  d->job_type_id = gauss_job_type_id(sym);
  /* system parameters */
  gauss_fchk_read_item("Number of atoms",TYPE_UINT,&d->n_atoms,f,file);
  gauss_fchk_read_item("Charge",TYPE_INT,&d->charge,f,file);
  gauss_fchk_read_item("Multiplicity",TYPE_UINT,&d->multiplicity,f,file);
  gauss_fchk_read_item("Number of electrons",TYPE_UINT,&d->n_electrons,f,file);
  gauss_fchk_read_item("Number of alpha electrons",TYPE_UINT,
    &d->n_alpha_electrons,f,file);
  gauss_fchk_read_item("Number of beta electrons",TYPE_UINT,
    &d->n_beta_electrons,f,file);
  /* atomic data */
  gauss_fchk_read_atom(d,f,file);
  /* basis functions */
  gauss_fchk_read_basis(d,f,file);
  /* molecular orbitals */
  gauss_fchk_read_mo(d,f,file);
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
