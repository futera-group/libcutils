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
#include <mol/molec.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Read structure from open XYZ file
  
   d - gromacs data structure
   f - open file stream */
void gmx_dat_fread_xyz(struct gmx_dat *d, FILE *f) {
  unsigned i,j;
  double mn,mx,dx;
  struct xyz_mol *x;
  /* input coordinates */
  x = mol_xyz_new();
  mol_xyz_fread(x,1,f);
  /* title */
  d->title = x->title;
  str_trim(d->title);
  x->title = NULL;
  /* compatibility with topology */
  if (d->n_frags && d->n_atoms!=x->n_atoms)
    msg_error_f("number of atoms in structure file inconsistent"
      " with topology (%d/%d)",1,x->n_atoms,d->n_atoms);
  /* save coordinates */
  if (d->n_frags) {
    d->crd = mat_ffree(d->crd,d->n_atoms);
    d->crd = mat_falloc(d->n_atoms,3);
    for (i=0; i<d->n_atoms; i++)
      for (j=0; j<3; j++)
        d->crd[i][j] = x->atom[i].coord[j];
    }
  /* save all data */
  else {
    d->n_atoms = x->n_atoms;
    d->n_resids = 1;
    d->n_mols = 1;
    d->n_frags = 1;
    d->crd = mat_falloc(d->n_atoms,3);
    /* molecules */
    d->mol = gmx_mol_new(d->n_mols);
    d->mol->name = str_copy_new("Molecule");
    d->mol->n_atoms = x->n_atoms;
    d->mol->n_resids = 1;
    /* residues */
    d->mol->res = gmx_res_new(d->mol->n_resids);
    d->mol->res->name = str_copy_new("MOL");
    d->mol->res->id = 0;
    d->mol->res->n_atoms = d->n_atoms;
    /* atoms */
    d->mol->res->atom = gmx_atom_new(d->n_atoms);
    for (i=0; i<d->n_atoms; i++) {
      d->mol->res->atom[i].id = i;
      d->mol->res->atom[i].num = x->atom[i].num;
      d->mol->res->atom[i].name = str_copy_new(x->atom[i].name);
      /* coordinates */
      for (j=0; j<3; j++)
        d->crd[i][j] = x->atom[i].coord[j];
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

/* Read structure from XYZ file

   d    - gromacs data structure
   name - name of the file */
void gmx_dat_read_xyz(struct gmx_dat *d, char *name) {
  FILE *f;
  f = file_open(name,"r");
  gmx_dat_fread_xyz(d,f);
  file_close(f);
  }

/* -------------------------------------------------------------------------- */
