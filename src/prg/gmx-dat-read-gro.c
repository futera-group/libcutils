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
#include <cmn/vector.h>
#include <mol/atom.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Read structure from open GRO file
  
   d - gromacs data structure
   f - open file stream */
void gmx_dat_fread_gro(struct gmx_dat *d, FILE *f) {
  unsigned i,j,k,l,ia,im,r0,r1,dm,n_atoms;
  char *line,res[80],sym[80],str[80];
  double crd[3],*x;
  struct gmx_atom *a;
  struct gmx_res *r = NULL;
  struct queue *q,*p,*t;
  /* title */
  d->title = str_read_line_new(f);
  if (!d->title)
    msg_error("structure file is empty",1);
  str_trim(d->title);
  /* number of atoms */
  line = str_read_line_new(f);
  if (sscanf(line,"%u",&n_atoms)!=1)
    msg_error_f("cannot read number of atoms in structure file",1);
  line = str_free(line);
  /* compatibility with topology */
  if (d->n_frags && d->n_atoms!=n_atoms)
    msg_error_f("number of atoms in structure file inconsistent"
      " with topology (%d/%d)",1,n_atoms,d->n_atoms);
  /* structure data */
  q = queue_alloc();
  t = queue_alloc();
  for (i=0; i<n_atoms; i++) {
    line = str_read_line_new(f);
    if (!line)
      msg_error("unexpected end of structure file",1);
    /* fixed-column format */
    if (sscanf(line,"%5u",&r1)!=1)
      msg_error_f("cannot read residuum ID of atom #%d in GRO file",1,i+1);
    if (sscanf(line+5,"%4s",res)!=1)
      msg_error_f("cannot read residuum name of atom #%d in GRO file",1,i+1);
    sprintf(str,"%6.6s",line+9);
    if (sscanf(str,"%6s",sym)!=1)
      msg_error_f("cannot read name of atom #%d in GRO file",1,i+1);
    if (sscanf(line+15,"%5u",&dm)!=1)
      msg_error_f("cannot read sequential number of atom #%d in GRO file",1,i+1);
    if (sscanf(line+20,"%8lf%8lf%8lf",&(crd[0]),&(crd[1]),&(crd[2]))!=3)
      msg_error_f("cannot read coordinates of atom #%d in GRO file",1,i+1);
    line = str_free(line);
    /* atoms */
    a = gmx_atom_new(1);
    a->num = atom_num_pdb(sym);
    a->name = str_copy_new(sym);
    /* residuum */
    if (!r || r0!=r1) {
      r = gmx_res_new(1);
      r->name = str_copy_new(res);
      p = queue_alloc();
      r->atom = (struct gmx_atom*)p;
      queue_add(q,r);
      }
    queue_add(t,vec_fcopy_new(crd,3));
    queue_add(p,a);
    r0 = r1;
    }
  /* save coordinates */
  if (d->n_frags) {
    d->crd = mat_ffree(d->crd,d->n_atoms);
    d->crd = mat_falloc(d->n_atoms,3);
    for (i=0,ia=0,im=0; i<d->n_frags; i++) {
      for (j=0; j<d->frag[i].n_rep; j++,im++) {
        for (k=0; k<d->frag[i].mol->n_resids; k++) {
          /* residuum */
          r = queue_get(q);
          if (!str_compare(d->frag[i].mol->res[k].name,r->name))
            msg_warn_f("different name of residuum #%d in molecule #%d"
              " (\"%s\" / \"%s\")",k+1,im+1,r->name,
              d->frag[i].mol->res[k].name);
          p = (struct queue*)r->atom;
          r->atom = NULL;
          for (l=0; l<d->frag[i].mol->res[k].n_atoms; l++,ia++) {
            /* atom type */
            a = queue_get(p);
            if (!str_compare(d->frag[i].mol->res[k].atom[l].name,a->name))
              msg_warn_f("different name of atom #%d in residuum #%d"
                " in molecule #%d (\"%s\" / \"%s\")",l+1,k+1,im+1,a->name,
                d->frag[i].mol->res[k].atom[l].name);
            gmx_atom_free(a,1);
            /* coordinate */
            x = queue_get(t);
            vec_fcopy_scaled(d->crd[ia],x,10.0,3);
            vec_ffree(x);
            }
          gmx_res_free(r,1);
          queue_free(p);
          }
        }
      }
    }
  /* save all data */
  else {
    d->n_atoms = n_atoms;
    d->n_resids = q->num;
    d->n_mols = 1;
    d->n_frags = 1;
    d->crd = mat_falloc(d->n_atoms,3);
    /* molecules */
    d->mol = gmx_mol_new(d->n_mols);
    d->mol->name = str_copy_new("Molecule");
    d->mol->n_atoms = n_atoms;
    d->mol->n_resids = q->num;
    /* residues */
    d->mol->res = gmx_res_new(d->mol->n_resids);
    for (i=0,ia=0; i<d->mol->n_resids; i++) {
      r = queue_get(q);
      p = (struct queue*)r->atom;
      r->atom = NULL;
      d->mol->res[i].name = str_copy_new(r->name);
      d->mol->res[i].id = i;
      d->mol->res[i].n_atoms = p->num;
      /* atoms */
      d->mol->res[i].atom = gmx_atom_new(d->mol->res[i].n_atoms);
      for (j=0; j<d->mol->res[i].n_atoms; j++,ia++) {
        a = queue_get(p);
        d->mol->res[i].atom[j].id = ia;
        d->mol->res[i].atom[j].num = atom_num_pdb(a->name);
        d->mol->res[i].atom[j].name = str_copy_new(a->name);
        d->mol->res[i].atom[j].mass = atom_mass(d->mol->res[i].atom[j].num);
        gmx_atom_free(a,1);
        /* coordinates */
        x = queue_get(t);
        vec_fcopy_scaled(d->crd[ia],x,10.0,3);
        vec_ffree(x);
        }
      gmx_res_free(r,1);
      queue_free(p);
      }
    /* compounds */
    d->frag = gmx_frag_new(d->n_frags);
    d->frag->name = str_copy_new("Molecule");
    d->frag->n_rep = 1;
    d->frag->mol = d->mol;
    }
  queue_free(q);
  queue_free(t);
  /* simulation box */
  line = str_read_line_new(f);
  if (line && sscanf(line,"%lf%lf%lf",&(crd[0]),&(crd[1]),&(crd[2]))==3) {
    d->box = vec_ffree(d->box);
    d->box = vec_fcopy_new(crd,3);
    vec_fscale(d->box,10.0,3);
    }
  }

/* Read structure from GRO file

   d    - gromacs data structure
   name - name of the file */
void gmx_dat_read_gro(struct gmx_dat *d, char *name) {
  FILE *f;
  f = file_open(name,"r");
  gmx_dat_fread_gro(d,f);
  file_close(f);
  }

/* -------------------------------------------------------------------------- */
