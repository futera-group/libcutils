/******************************************************************************\
 *                                                                            * 
 *  Libcutils - library of C function                                         * 
 *                                                                            *
 *  Version:             3.4                                                  * 
 *  Date:                26/03/2020                                           *
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
#include <cmn/string.h>
#include <cmn/vector.h>
#include <mol/atom.h>
#include <mol/molec.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Read structure from open PDB file
  
   d - gromacs data structure
   f - open file stream */
void gmx_dat_fread_pdb(struct gmx_dat *d, FILE *f) {
  unsigned i,j,k,ia,ir,n_atoms;
  double dx,mn,mx;
  struct pdb_mol *p;
  /* input coordinates */
  p = mol_pdb_new();
  mol_pdb_fread(p,1,f);
  /* title */
  d->title = p->title;
  str_trim(d->title);
  p->title = NULL;
  /* compatibility with topology */
  n_atoms = mol_pdb_atom_n(p);
  if (d->n_frags) {
    if (d->n_atoms!=n_atoms)
      msg_error_f("number of atoms in structure file inconsistent"
        " with topology (%d/%d)",1,n_atoms,d->n_atoms);
    if (d->n_resids!=p->n_res) 
      msg_error_f("number of residues in structure file inconsistent"
        " with topology (%d/%d)",1,p->n_res,d->n_resids);
    ir = 0;
    for (i=0; i<d->n_frags; i++) {
      for (j=0; j<d->frag[i].n_rep; j++) {
        for (k=0; k<d->frag[i].mol->n_resids; k++) {
          if (d->frag[i].mol->res[k].n_atoms != p->res[ir].n_atoms)
            msg_error_f("number of atoms in residuum #%d is inconsistent"
              " with topology (%d/%d)",1,ir+1,p->res[ir].n_atoms,
              d->frag[i].mol->res[k].n_atoms);
          ir++;
          }
        }
      }
    }
  /* save coordinates */
  if (d->n_frags) {
    d->crd = mat_ffree(d->crd,d->n_atoms);
    d->crd = mat_falloc(d->n_atoms,3);
    for (i=0,ia=0; i<p->n_res; i++)
      for (j=0; j<p->res[i].n_atoms; j++,ia++)
        for (k=0; k<3; k++)
          d->crd[ia][k] = p->res[i].atom[j].coord[k];
    }
  /* save all data */
  else {
    d->n_atoms = n_atoms;
    d->n_resids = p->n_res;
    d->n_mols = 1;
    d->n_frags = 1;
    d->crd = mat_falloc(d->n_atoms,3);
    /* molecules */
    d->mol = gmx_mol_new(d->n_mols);
    d->mol->name = str_copy_new("Molecule");
    d->mol->n_atoms = n_atoms;
    d->mol->n_resids = p->n_res;
    /* residues */
    d->mol->res = gmx_res_new(d->mol->n_resids);
    for (i=0,ia=0; i<p->n_res; i++) {
      d->mol->res[i].name = str_copy_new(p->res[i].name);
      d->mol->res[i].id = p->res[i].id;
      d->mol->res[i].n_atoms = p->res[i].n_atoms;
      /* atoms */
      d->mol->res[i].atom = gmx_atom_new(d->mol->res[i].n_atoms);
      for (j=0; j<d->mol->res[i].n_atoms; j++,ia++) {
        d->mol->res[i].atom[j].id = p->res[i].atom[j].id;
        d->mol->res[i].atom[j].num = atom_num(p->res[i].atom[j].name);
        d->mol->res[i].atom[j].name = str_copy_new(p->res[i].atom[j].name);
        /* coordinates */
        for (k=0; k<3; k++)
          d->crd[ia][k] = p->res[i].atom[j].coord[k];
        }
      }
    /* compounds */
    d->frag = gmx_frag_new(d->n_frags);
    d->frag->name = str_copy_new("Molecule");
    d->frag->n_rep = 1;
    d->frag->mol = d->mol;
    }
  /* simulation box */
  if (d->n_atoms) {
    d->box = vec_ffree(d->box);
    d->box = vec_falloc(3);
    dx = 5.0;
    for (i=0; i<3; i++) {
      mn =  9.9e+99;
      mx = -9.9e+99;
      for (j=0; j<d->n_atoms; j++) {
        if (d->crd[j][i] < mn)
          mn = d->crd[j][i];
        if (d->crd[j][i] > mx)
          mx = d->crd[j][i];
        d->box[i] = mx-mn+dx;
        }
      }
    }
  }

/* Read structure from PDB file

   d    - gromacs data structure
   name - name of the file */
void gmx_dat_read_pdb(struct gmx_dat *d, char *name) {
  FILE *f;
  f = file_open(name,"r");
  gmx_dat_fread_pdb(d,f);
  file_close(f);
  }

/* -------------------------------------------------------------------------- */
