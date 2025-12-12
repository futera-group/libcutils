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

#include <math.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "mol/aacid.h"
#include "mol/atom.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* Shift structure geometrical center to the coordinate origin
 
   p - the structure data */
void aacid_chain_center(struct pdb_mol *p) {
  unsigned i,j,k,n = 0;
  double c[3];
  vec_fset(c,0.0,3);
  for (i=0; i<p->n_res; i++)
    for (j=0; j<p->res[i].n_atoms; j++) {
      for (k=0; k<3; k++)
        c[k] += p->res[i].atom[j].coord[k];
      n++;
      }
  if (n) {
    for (i=0; i<3; i++)
      c[i] /= n;
    for (i=0; i<p->n_res; i++)
      for (j=0; j<p->res[i].n_atoms; j++)
        for (k=0; k<3; k++)
          p->res[i].atom[j].coord[k] -= c[k];
    }
  }

/* -------------------------------------------------------------------------- */

/* Delete all hydrogen atom from the structure
 
   p - the structure data */
void aacid_chain_del_h(struct pdb_mol *p) {
  unsigned i,j,k,n;
  struct pdb_atom *a;
  /* delete hydrogens */
  for (i=0; i<p->n_res; i++) {
    /* count heavy atoms */
    for (j=0,n=0; j<p->res[i].n_atoms; j++)
      if (atom_num_pdb(p->res[i].atom[j].name)!=1)
        n++;
    /* create new atom array */
    if (n<p->res[i].n_atoms) {
      a = mol_pdb_atom_new(n);
      for (j=0,k=0; j<p->res[i].n_atoms; j++)
        if (atom_num_pdb(p->res[i].atom[j].name)!=1) {
          a[k].name = str_copy_new(p->res[i].atom[j].name);
          vec_fcopy(a[k].coord,p->res[i].atom[j].coord,3);
          k++;
          }
      /* replace the atom array */
      mol_pdb_atom_free(p->res[i].atom,p->res[i].n_atoms);
      p->res[i].n_atoms = n;
      p->res[i].atom = a;
      }
    }
  /* update number of atom IDs */
  for (i=0,k=0; i<p->n_res; i++)
    for (j=0; j<p->res[i].n_atoms; j++)
      p->res[i].atom[j].id = k++;
  }

/* -------------------------------------------------------------------------- */

/* Add amino acid to peptide chain 
 
   d  - molecular data structure (PDB)
   id - sequential ID
   z  - amino-acid structure (Z-matrix) */
void aacid_chain_add(struct pdb_mol *d, unsigned id, struct list *z) {
  unsigned i,i2,i3,i4;
  double v,r[3];
  struct ldata *p;
  struct pdb_res *t0,*t1;
  struct zmt_atom *a;
  /* amino acid data */
  t0 = (id ? d->res+id-1 : NULL);
  t1 = d->res+id;
  /* save atomic data */
  t1->n_atoms = z->num;
  t1->atom = mol_pdb_atom_new(t1->n_atoms);
  for (p=z->first,i=0; p; p=p->l_next,i++) {
    a = (struct zmt_atom*)p->l_data;
    /* atomic name */
    t1->atom[i].name = str_copy_new(a->name);
    /* coordinates */
    vec_fset(r,0.0,3);
    /* first atom */
    if (i==0) {
      if (id) {
        i2 = mol_pdb_res_atom_find(t0," C  ",NULL);
        i3 = mol_pdb_res_atom_find(t0," O  ",NULL);
        i4 = mol_pdb_res_atom_find(t0," CA ",NULL);
        mol_zmt_atom_vec(r,t0->atom[i2].coord,t0->atom[i3].coord,
          t0->atom[i4].coord,1.32,120.0,180.0);
        vec_fadd(t1->atom[i].coord,t0->atom[i2].coord,r,3);
        }
      else
        vec_fset(t1->atom[i].coord,0.0,3);
      }
    /* second atom */
    else if (i==1) {
      if (id) {
        i3 = mol_pdb_res_atom_find(t0," C  ",NULL);
        i4 = mol_pdb_res_atom_find(t0," O  ",NULL);
        mol_zmt_atom_vec(r,t1->atom[a->bond_id].coord,
          t0->atom[i3].coord,t0->atom[i4].coord,a->bond_val,120.0,180.0);
        }
      else
        r[0] = a->bond_val;
      vec_fadd(t1->atom[i].coord,t1->atom[a->bond_id].coord,r,3);
      }
    /* third atom */
    else if (i==2) {
      if (id) {
        i4 = mol_pdb_res_atom_find(t0," O  ",NULL);
        mol_zmt_atom_vec(r,t1->atom[a->bond_id].coord,
          t1->atom[a->angle_id].coord,t0->atom[i4].coord,
          a->bond_val,a->angle_val,t1->id==AA_PRO ? 0.0 : 180.0);
        } 
      else {
        v = (a->angle_val<90.0 ? 180.0-a->angle_val : a->angle_val);
        r[0] = a->bond_val*cos(M_PI*v/180.0);
        r[1] = a->bond_val*sin(M_PI*v/180.0);
        }
      vec_fadd(t1->atom[i].coord,t1->atom[a->bond_id].coord,r,3);
      }
    /* other atoms */
    else {
      mol_zmt_atom_vec(r,
        t1->atom[a->bond_id].coord,
        t1->atom[a->angle_id].coord,
        t1->atom[a->dihed_id].coord, 
        a->bond_val,a->angle_val,a->dihed_val);
      vec_fadd(t1->atom[i].coord,t1->atom[a->bond_id].coord,r,3);
      }
    }
  }

/* create chain of amino-acids (peptide structure)
 
   seq  - the amino-acid sequence (one-letter notation)
   form - protonation form of amino acids
   hydr - presence of hydrogens */
struct pdb_mol* aacid_chain(char *seq, short form, short hydr) {
  unsigned i,j,n,id;
  struct list *t;
  struct pdb_mol *p;
  /* prepare data structure */
  p = mol_pdb_new();
  p->n_res = str_length(seq);
  p->res = mol_pdb_res_new(p->n_res);
  /* convert symbolic sequence to amino-acid data array */
  for (i=0,n=0; i<p->n_res; i++) {
    p->res[i].id = i+1;
    id = aacid_id_1(seq[i]);
    p->res[i].name = str_copy_new(str_upcase(aacid_name_3(id)));
    /* single amino acid */
    if (p->n_res==1)
      t = aacid_geom_def(id,form);
    /* N-terminus */
    else if (i==0)
      t = aacid_geom_def(id,form==AA_FORM_NEUTRAL ? 
            AA_FORM_NTER_N : AA_FORM_NTER_P);
    /* C-terminus */
    else if ((i+1)==p->n_res)
      t = aacid_geom_def(id,form==AA_FORM_NEUTRAL ? 
            AA_FORM_CTER_N : AA_FORM_CTER_P);
    /* chain amino acid */
    else
      t = aacid_geom_def(id,form==AA_FORM_NEUTRAL ? 
            AA_FORM_CHAIN_N : AA_FORM_CHAIN_P);
    /* add amino acid to peptide chain */
    aacid_chain_add(p,i,t);
    for (j=0; j<p->res[i].n_atoms; j++)
      p->res[i].atom[j].id = n++;
    }
  /* delete hydrogens */
  if (!hydr)
    aacid_chain_del_h(p);
  /* shift structure to geometrical center */
  aacid_chain_center(p);
  return(p);
  }

/* -------------------------------------------------------------------------- */
