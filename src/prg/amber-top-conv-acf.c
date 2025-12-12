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

/* convert amber topology and coordinates to acf molecular file format
 
   t - pointer to amber topology struct
   c - pointer to array with coordinates */
struct acf_mol *amber_top_acf(struct amber_top *t, double *c) {
  unsigned i,j,k,*r,nr,a0 = 0,a1 = 0,id = 0;
  double chrg = 0.0;
  struct acf_mol *m;
  m = mol_acf_new();
  r = amber_res_unique(t,&nr);
  for (i=0; i<nr; i++)
    m->n_atoms += amber_res_natoms(t,r[i]);
  m->atom = mol_acf_atom_new(m->n_atoms);
  for (i=0; i<nr; i++) {
    amber_res_atom_id(t,r[i],&a0,&a1);
    for (j=a0; j<a1; j++) {
      m->atom[id].name = str_copy_new(t->atom_names[j]);
      m->atom[id].type = str_copy_new(t->atom_types[j]);
      m->atom[id].num = atom_num_pdb(t->atom_names[j]);
      m->atom[id].charge = t->charge[j]/AMBER_TOP_CHARGE;
      chrg += m->atom[id].charge;
      for (k=0; k<3; k++)
        m->atom[id].coord[k] = c[3*j+k];
      id++;
      }
    }
  m->charge = chrg;
  m->name = str_copy_new("MOL");
  mol_acf_set_formula(m);
  mol_acf_set_bonds(m);
  return(m);
  }

/* convert amber topology and coordinates to acf molecular file format
 
   t - pointer to amber topology struct
   c - pointer to array with coordinates */
struct acf_mol *amber_top_acf_all(struct amber_top *t, double *c) {
  double chrg = 0.0;
  unsigned i,j;
  struct acf_mol *m;
  m = mol_acf_new();
  m->n_atoms = t->pointers[AMBER_POINTER_NATOM];
  m->atom = mol_acf_atom_new(m->n_atoms);
  for (i=0; i<m->n_atoms; i++) {
    m->atom[i].name = str_copy_new(t->atom_names[i]);
    m->atom[i].type = str_copy_new(t->atom_types[i]);
    m->atom[i].num = atom_num_pdb(t->atom_names[i]);
    m->atom[i].charge = t->charge[i]/AMBER_TOP_CHARGE;
    chrg += m->atom[i].charge;
    for (j=0; j<3; j++)
      m->atom[i].coord[j] = c[3*i+j];
    }
  m->charge = chrg;
  m->name = str_copy_new("MOL");
  mol_acf_set_formula(m);
  mol_acf_set_bonds(m);
  return(m);
  }
/* -------------------------------------------------------------------------- */
