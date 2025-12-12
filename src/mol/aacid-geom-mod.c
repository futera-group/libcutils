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

#include <stdlib.h>
#include <cmn/message.h>
#include <cmn/vector.h>
#include "mol/aacid.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* Delete extra atoms from amino-acid structure

   x    - structure of single amino acid (XYZ)
   form - structure form specification */
void aacid_geom_mod_del(struct xyz_mol *x, short form) {
  switch (form) {
    /* neutral N-terminus */
    case AA_FORM_NTER_N:
      mol_xyz_atom_del_name(x," OH ");
      mol_xyz_atom_del_name(x," HO ");
      break;
    /* neutral C-terminus */
    case AA_FORM_CTER_N: 
      mol_xyz_atom_del_name(x," HN ");
      break;
    /* neutral peptide chain */
    case AA_FORM_CHAIN_N:
      mol_xyz_atom_del_name(x," HN ");
      mol_xyz_atom_del_name(x," OH ");
      mol_xyz_atom_del_name(x," HO ");
      break;
    /* protonated N-terminus */
    case AA_FORM_NTER_P:
      mol_xyz_atom_del_name(x," OXT");
      break;
    /* deprotonated C-terminus */
    case AA_FORM_CTER_P:
      mol_xyz_atom_del_name(x," HN1");
      mol_xyz_atom_del_name(x," HN2");
      break;
    /* de/protonated peptide chain */
    case AA_FORM_CHAIN_P:
      mol_xyz_atom_del_name(x," HN1");
      mol_xyz_atom_del_name(x," HN2");
      mol_xyz_atom_del_name(x," OXT");
      break;
    }
  }

/* -------------------------------------------------------------------------- */

/* Set atom ID for reordering array 
 
   x    - the amino-acid structure
   aa   - amino acid ID
   name - name of the atom */
unsigned aacid_geom_mod_order_id(struct xyz_mol *x, short aa, char* name) {
  unsigned id;
  short found;
  id = mol_xyz_atom_find(x,name,&found);
  if (!found)
    msg_error_f("cannot file atom \"%s\" in %s\n",
      1,name,aacid_name_full(aa));
  return(id);
  }

/* Reorder atoms in amino-acid structure 
 
   x    - the amino-acid structure
   aa   - amino acid ID
   form - structure form (N-terminus,C-terminus,chain) */
void aacid_geom_mod_order(struct xyz_mol *x, short aa, short form) {
  short found;
  unsigned i,j,n,*r;
  /* atom ID array */
  n = x->n_atoms;
  r = vec_ualloc(n);
  vec_uset(r,n+1,n);
  switch (form) {
    /* neutral N-terminus */
    case AA_FORM_NTER_N:
      if (aa==AA_PRO) {
        r[0] = aacid_geom_mod_order_id(x,aa," N  ");
        r[1] = aacid_geom_mod_order_id(x,aa," HN ");
        }
      else {
        r[0] = aacid_geom_mod_order_id(x,aa," N  ");
        r[1] = aacid_geom_mod_order_id(x,aa," H  ");
        r[2] = aacid_geom_mod_order_id(x,aa," HN ");
        }
      r[n-2] = aacid_geom_mod_order_id(x,aa," C  ");
      r[n-1] = aacid_geom_mod_order_id(x,aa," O  ");
      break;
    /* neutral C-terminus */
    case AA_FORM_CTER_N: 
      if (aa==AA_PRO)
        r[0] = aacid_geom_mod_order_id(x,aa," N  ");
      else {
        r[0] = aacid_geom_mod_order_id(x,aa," N  ");
        r[1] = aacid_geom_mod_order_id(x,aa," H  ");
        }
      r[n-4] = aacid_geom_mod_order_id(x,aa," C  ");
      r[n-3] = aacid_geom_mod_order_id(x,aa," O  ");
      r[n-2] = aacid_geom_mod_order_id(x,aa," OH ");
      r[n-1] = aacid_geom_mod_order_id(x,aa," HO ");
      break;
    /* protonated N-terminus */
    case AA_FORM_NTER_P:
      if (aa==AA_PRO) {
        r[0] = aacid_geom_mod_order_id(x,aa," N  ");
        r[1] = aacid_geom_mod_order_id(x,aa," HN1");
        r[2] = aacid_geom_mod_order_id(x,aa," HN2");
        }
      else {
        r[0] = aacid_geom_mod_order_id(x,aa," N  ");
        r[1] = aacid_geom_mod_order_id(x,aa," H  ");
        r[2] = aacid_geom_mod_order_id(x,aa," HN1");
        r[3] = aacid_geom_mod_order_id(x,aa," HN2");
        }
      r[n-2] = aacid_geom_mod_order_id(x,aa," C  ");
      r[n-1] = aacid_geom_mod_order_id(x,aa," O  ");
      break;
    /* deprotonated C-terminus */
    case AA_FORM_CTER_P:
      if (aa==AA_PRO)
        r[0] = aacid_geom_mod_order_id(x,aa," N  ");
      else {
        r[0] = aacid_geom_mod_order_id(x,aa," N  ");
        r[1] = aacid_geom_mod_order_id(x,aa," H  ");
        }
      r[n-3] = aacid_geom_mod_order_id(x,aa," C  ");
      r[n-2] = aacid_geom_mod_order_id(x,aa," O  ");
      r[n-1] = aacid_geom_mod_order_id(x,aa," OXT");
      break;
    /* peptide chain */
    case AA_FORM_CHAIN_N:
    case AA_FORM_CHAIN_P:
      if (aa==AA_PRO)
        r[0] = aacid_geom_mod_order_id(x,aa," N  ");
      else {
        r[0] = aacid_geom_mod_order_id(x,aa," N  ");
        r[1] = aacid_geom_mod_order_id(x,aa," H  ");
        }
      r[n-2] = aacid_geom_mod_order_id(x,aa," C  ");
      r[n-1] = aacid_geom_mod_order_id(x,aa," O  ");
      break;
    }
  /* reordering IDs */
  for (i=0; i<n; i++) {
    for (j=0,found=0; j<n; j++)
      if (r[j]==i) {
        found = 1;
        break;
        }
    if (!found) {
      for (j=0; j<n; j++)
        if (r[j]>n) {
          r[j] = i;
          break;
          }
      }
    }
  mol_xyz_atom_order(x,r);
  vec_ufree(r);
  }

/* -------------------------------------------------------------------------- */

/* Structure modifications (C-terminus,N-terminus,chain)

   z    - structure of single amino acid (Z-matrix)
   aa   - amino acid ID
   form - structure form specification */
void aacid_geom_mod(struct list *z, short aa, short form) {
  unsigned i;
  struct ldata *p;
  struct zmt_mol *t;
  struct zmt_atom *a;
  struct xyz_mol *x;
  /* create Z-matrix */
  t = mol_zmt_new();
  t->n_atoms = z->num;
  t->atom = mol_zmt_atom_new(t->n_atoms);
  for (p=z->first,i=0; p; p=p->l_next,i++)
    mol_zmt_atom_copy((struct zmt_atom*)p->l_data,t->atom+i);
  /* conversion to XYZ */
  x = mol_zmt_xyz(t);
  mol_zmt_free(t);
  /* delete unwanted atoms */
  aacid_geom_mod_del(x,form);
  /* atom reorder */
  aacid_geom_mod_order(x,aa,form);
  /* conversion to ZMT */
  t = mol_xyz_zmt(x,NULL);
  mol_xyz_free(x);
  /* update Z-matrix */
  list_clean(z,free);
  for (i=0; i<t->n_atoms; i++) {
    a = mol_zmt_atom_new(1);
    mol_zmt_atom_copy(t->atom+i,a);
    list_add_end_p(z,a);
    }
  mol_zmt_free(t);
  }

/* -------------------------------------------------------------------------- */
