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
#include <cmn/string.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* allocate new xyz molecule struct */
struct xyz_mol *mol_xyz_new(void) {
  struct xyz_mol *x = NULL;
  x = (struct xyz_mol*)malloc(sizeof(struct xyz_mol));
  if (!x)
    msg_error("cannot allocate memory for xyz molecular data struct",1);
  /* initialization */
  x->title = NULL;
  x->n_atoms = 0;
  x->atom = NULL;
  return(x);
  }

/* allocate array of xyz atom data struct
 
   n - number of atoms */
struct xyz_atom *mol_xyz_atom_new(unsigned n) {
  unsigned i,j;
  struct xyz_atom *a = NULL;
  /* allocate memory */
  a = (struct xyz_atom*)malloc(n*sizeof(struct xyz_atom));
  if (!a)
    msg_error("cannot allocate memory of xyz atom data array",1);
  /* initialization */
  for (i=0; i<n; i++) {
    a[i].name = NULL;
    a[i].num = 0;
    for (j=0; j<3; j++)
      a[i].coord[j] = 0.0;
    }
  return(a);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for xyz molecule struct

   x - pointer to xyz molecular data struct */
void mol_xyz_free(struct xyz_mol *x) {
  if (x) {
    mol_xyz_atom_free(x->atom,x->n_atoms);
    free(x);
    }
  }

/* free memory allocated for xyz atom array

   a - the atom array
   n - number of atoms */
void mol_xyz_atom_free(struct xyz_atom *a, unsigned n) {
  unsigned i;
  if (a) {
    for (i=0; i<n; i++)
      a[i].name = str_free(a[i].name);
    free(a);
    }
  }

/* -------------------------------------------------------------------------- */

/* copy XYZ molecular data struct
  
   x - pointer to XYZ molecular data struct */
struct xyz_mol *mol_xyz_copy_new(struct xyz_mol *x) {
  unsigned i,j;
  struct xyz_mol *y;
  /* new xyz structure */
  y = mol_xyz_new();
  y->n_atoms = x->n_atoms;
  y->atom = mol_xyz_atom_new(y->n_atoms);
  y->title = str_copy_new(x->title);
  for (i=0; i<x->n_atoms; i++) {
    y->atom[i].name = str_copy_new(x->atom[i].name);
    y->atom[i].num = x->atom[i].num;
    for (j=0; j<3; j++)
      y->atom[i].coord[j] = x->atom[i].coord[j];
    }
  return(y);
  }

/* -------------------------------------------------------------------------- */
