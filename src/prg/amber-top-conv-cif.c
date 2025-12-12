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
#include <mol/cell.h>
#include <mol/molec.h>
#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

/* set box parameters for cif structure converted from amber topology
 
   t - pointer to amber topology struct
   c - pointer to cif molecular data struct */
void amber_top_cif_box(struct amber_top *t, struct cif_mol *c) {
  double max[3] = {-9.0E+90,-9.0E+90,-9.0E+90};
  double min[3] = { 9.0E+90, 9.0E+90, 9.0E+90};
  unsigned i,j;
  if (t->pointers[AMBER_POINTER_BOX]) 
    cell_set_side_angle_v(c->cell,t->box_params+3,t->box_params);
  else {
    for (i=0; i<c->n_atoms; i++)
      for (j=0; j<3; j++) {
        if (c->atom[i].coord[j]<min[j])
          min[j] = c->atom[i].coord[j];
        if (c->atom[i].coord[j]>max[j])
          max[j] = c->atom[i].coord[j];
        }
    cell_set_side_angle(c->cell,
      max[0]-min[0],max[1]-min[1],max[2]-min[2],
      90.0,90.0,90.0);
    }
  c->cell->space_group_name = str_copy_new("P 1");
  c->cell->space_group_id = 1;
  }

/* convert amber topology and coordinates to cif molecular file format
 
   t - pointer to amber topology struct
   x - pointer to array with coordinates */
struct cif_mol *amber_top_cif(struct amber_top *t, double *x) {
  unsigned i,j;
  struct cif_mol *c;
  /* cif structure */
  c = mol_cif_new();
  c->n_atoms = t->pointers[AMBER_POINTER_NATOM];
  c->atom = mol_cif_atom_new(c->n_atoms);
  for (i=0; i<c->n_atoms; i++) {
    for (j=0; j<3; j++)
      c->atom[i].coord[j] = x[3*i+j];
    c->atom[i].num = atom_num_pdb(t->atom_names[i]);
    c->atom[i].name = str_copy_new(t->atom_names[i]);
    str_trim(c->atom[i].name);
    }
  amber_top_cif_box(t,c);
  return(c);
  }

/* convert amber topology and coordinates to cif molecular file format
 
   t - pointer to amber topology struct
   x - pointer to array with coordinates
   r - residuum ID (0 for the whole structure) */
struct cif_mol *amber_top_cif_r(struct amber_top *t, double *x, unsigned r) {
  unsigned i,j,a0 = 0,a1 = 0;
  struct cif_mol *c = NULL;
  if (!r)
    c = amber_top_cif(t,x);
  else {
    c = mol_cif_new();
    amber_res_atom_id(t,r-1,&a0,&a1);
    c->n_atoms = a1-a0;
    c->atom = mol_cif_atom_new(c->n_atoms);
    for (i=0; i<c->n_atoms; i++) {
      for (j=0; j<3; j++)
        c->atom[i].coord[j] = x[3*(a0+i)+j];
      c->atom[i].num = atom_num_pdb(t->atom_names[a0+i]);
      c->atom[i].name = str_copy_new(t->atom_names[a0+i]);
      str_trim(c->atom[i].name);
      }
    amber_top_cif_box(t,c);
    }
  return(c);
  }

/* convert amber topology and coordinates to cif molecular file format
 
   t - pointer to amber topology struct
   x - pointer to array with coordinates
   v - array with residuum IDs
   n - number of residuum IDs */
struct cif_mol *amber_top_cif_rv(struct amber_top *t, double *x,
  unsigned *v, unsigned n) {
  unsigned i,j,k,nt = 0,id = 0,a0 = 0,a1 = 0;
  struct cif_mol *c = NULL;
  if (!v || !n)
    c = amber_top_cif(t,x);
  else {
    /* total number of atom */
    for (i=0; i<n; i++)
      nt += amber_res_natoms(t,v[i]-1);
    /* cif structure */
    c = mol_cif_new();
    c->n_atoms = nt;
    c->atom = mol_cif_atom_new(c->n_atoms);
    for (i=0; i<n; i++) {
      amber_res_atom_id(t,v[i]-1,&a0,&a1);
      for (j=a0; j<a1; j++) {
        for (k=0; k<3; k++)
          c->atom[id].coord[k] = x[3*j+k];
        c->atom[id].num = atom_num_pdb(t->atom_names[j]);
        c->atom[id].name = str_copy_new(t->atom_names[j]);
        str_trim(c->atom[id].name);
        id++;
        }
      }
    amber_top_cif_box(t,c);
    }
  return(c);
  }

/* -------------------------------------------------------------------------- */
