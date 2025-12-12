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

#include <cmn/string.h>
#include <mol/atom.h>
#include <mol/molec.h>
#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

/* convert amber topology and coordinates to xyz molecular file format
 
   t - pointer to amber topology struct
   c - pointer to array with coordinates */
struct xyz_mol *amber_top_xyz(struct amber_top *t, double *c) {
  unsigned i,j;
  struct xyz_mol *x;
  /* xyz structure */
  x = mol_xyz_new();
  x->title = str_copy_new("The whole structure");
  x->n_atoms = t->pointers[AMBER_POINTER_NATOM];
  x->atom = mol_xyz_atom_new(x->n_atoms);
  for (i=0; i<x->n_atoms; i++) {
    for (j=0; j<3; j++)
      x->atom[i].coord[j] = c[3*i+j];
    x->atom[i].num = atom_num_pdb(t->atom_names[i]);
    }
  return(x);
  }

/* convert amber topology and coordinates to xyz molecular file format
 
   t - pointer to amber topology struct
   c - pointer to array with coordinates
   r - residuum ID (0 for the whole structure) */
struct xyz_mol *amber_top_xyz_r(struct amber_top *t, double *c, unsigned r) {
  char text[1024];
  unsigned i,j,id = 0,a0 = 0,a1 = 0;
  struct xyz_mol *x = NULL;
  if (!r)
    x = amber_top_xyz(t,c);
  else {
    x = mol_xyz_new();
    sprintf(text,"Residuum #%d: %s",r,t->res_names[r-1]);
    x->title = str_copy_new(text);
    amber_res_atom_id(t,r-1,&a0,&a1);
    x->n_atoms = a1-a0;
    x->atom = mol_xyz_atom_new(x->n_atoms);
    for (i=a0; i<a1; i++) {
      for (j=0; j<3; j++)
        x->atom[id].coord[j] = c[3*i+j];
      x->atom[id].num = atom_num_pdb(t->atom_names[i]);
      id++;
      }
    }
  return(x);
  }

/* convert amber topology and coordinates to xyz molecular file format
 
   t - pointer to amber topology struct
   c - pointer to array with coordinates
   v - array with residuum IDs
   n - number of residuum IDs */
struct xyz_mol *amber_top_xyz_rv(struct amber_top *t, double *c,
  unsigned *v, unsigned n) {
  unsigned i,j,k,nt = 0,id = 0,a0 = 0,a1 = 0;
  char range[1024],text[2048];
  struct xyz_mol *x = NULL;
  /* the whole structure */
  if (!v || !n)
    x = amber_top_xyz(t,c);
  else {
    str_unum_range(range,v,n);
    /* total number of atom */
    for (i=0; i<n; i++)
      nt += amber_res_natoms(t,v[i]-1);
    /* xyz structure */
    x = mol_xyz_new();
    x->n_atoms = nt;
    x->atom = mol_xyz_atom_new(x->n_atoms);
    sprintf(text,"Residui #%s",range);
    x->title = str_copy_new(text);
    /* atom conversion */
    for (i=0; i<n; i++) {
      amber_res_atom_id(t,v[i]-1,&a0,&a1);
      for (j=a0; j<a1; j++) {
        for (k=0; k<3; k++)
          x->atom[id].coord[k] = c[3*j+k];
        x->atom[id].num = atom_num_pdb(t->atom_names[j]);
        id++;
        }
      }
    }
  return(x);
  }

/* -------------------------------------------------------------------------- */
