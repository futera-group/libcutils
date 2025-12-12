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

#include <cmn/message.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* order atoms according to new indices

   x - pointer to XYZ molecular data struct
   v - vector with atomic indices */
void mol_xyz_atom_order(struct xyz_mol *x, unsigned *v) {
  unsigned i,j;
  struct xyz_mol *y;
  y = mol_xyz_copy_new(x);
  for (i=0; i<x->n_atoms; i++) {
    x->atom[i].name = str_copy_new(y->atom[v[i]].name);
    x->atom[i].num = y->atom[v[i]].num;
    for (j=0; j<3; j++)
      x->atom[i].coord[j] = y->atom[v[i]].coord[j];
    }
  mol_xyz_free(y);
  }

/* -------------------------------------------------------------------------- */

/* delete specified atom from the structure
 
   x  - pointer to XYZ molecular data struct
   id - sequential ID of the atom */
void mol_xyz_atom_del(struct xyz_mol *x, unsigned id) {
  unsigned i;
  struct xyz_mol *y;
  /* create new structure */
  y = mol_xyz_new();
  y->n_atoms = x->n_atoms-1;
  y->atom = mol_xyz_atom_new(y->n_atoms);
  /* omit the deleted atom */
  for (i=0; i<id; i++) {
    y->atom[i].name = str_copy_new(x->atom[i].name);
    y->atom[i].num = x->atom[i].num;
    vec_fcopy(y->atom[i].coord,x->atom[i].coord,3);
    }
  for (i=id+1; i<x->n_atoms; i++) {
    y->atom[i-1].name = str_copy_new(x->atom[i].name);
    y->atom[i-1].num = x->atom[i].num;
    vec_fcopy(y->atom[i-1].coord,x->atom[i].coord,3);
    }
  /* update the original structure */
  mol_xyz_atom_free(x->atom,x->n_atoms);
  x->n_atoms = y->n_atoms;
  x->atom = y->atom;
  /* clean memory */
  y->n_atoms = 0;
  y->atom = NULL;
  mol_xyz_free(y);
  }

/* Delete atom with specified name from the structure
 
   x    - the XYZ data structure
   name - name of the atom */
void mol_xyz_atom_del_name(struct xyz_mol *x, char *name) {
  unsigned id;
  short found;
  id = mol_xyz_atom_find(x,name,&found);
  if (!found)
    msg_error_f("atom \"%s\" not found in xyz structure",1,name);
  mol_xyz_atom_del(x,id);
  }

/* -------------------------------------------------------------------------- */

/* Find atom in the XYZ data structure and return its ID
 
   x  - the XYZ data structure
   a  - name of the atom
   ok - status (optional) */
unsigned mol_xyz_atom_find(struct xyz_mol *x, char *a, short *ok) {
  unsigned i,id = 0;
  if (ok)
    (*ok) = 0;
  for (i=0; i<x->n_atoms; i++) {
    if (str_compare(x->atom[i].name,a)) {
      if (ok)
        (*ok) = 1;
      id = i;
      break;
      }
    }
  return(id);
  }

/* -------------------------------------------------------------------------- */
