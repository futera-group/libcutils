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

/* allocate new cif molecule struct */
struct cif_mol* mol_cif_new(void) {
  struct cif_mol *c = NULL;
  c = (struct cif_mol*)malloc(sizeof(struct cif_mol));
  if (!c)
    msg_error("cannot allocate memory for cif molecular data struct",1);
  /* initialization */
  c->n_atoms = 0;
  c->cell = cell_new();
  c->atom = NULL;
  return(c);
  }

/* allocate array of cif atom data struct
 
   n - number of atoms */
struct cif_atom* mol_cif_atom_new(unsigned n) {
  unsigned i,j;
  struct cif_atom *a = NULL;
  a = (struct cif_atom*)malloc(n*sizeof(struct cif_atom));
  if (!a)
    msg_error("cannot allocate memory of cif atom data array",1);
  /* initialization */
  for (i=0; i<n; i++) {
    a[i].num = 0;
    a[i].name = NULL;
    for (j=0; j<3; j++)
      a[i].coord[j] = 0.0;
    }
  return(a);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for cif molecule struct
 
   c - pointer to cif molecular data struct */
void mol_cif_free(struct cif_mol *c) {
  if (c) {
    cell_free(c->cell);
    mol_cif_atom_free(c->atom,c->n_atoms);
    }
  }

/* free memory allocated for array of cif atom data structures 

   c - array of atoms
   n - number of atoms */
void mol_cif_atom_free(struct cif_atom *c, unsigned n) {
  unsigned i;
  if (c) {
    for (i=0; i<n; i++)
      str_free(c[i].name);
    free(c);
    }
  }

/* -------------------------------------------------------------------------- */

/* copy cif molecular data struct
 
   c - pointer to cif molecular data struct */
struct cif_mol* mol_cif_copy_new(struct cif_mol *c) {
  unsigned i;
  struct cif_mol *d;
  d = mol_cif_new();
  d->n_atoms = c->n_atoms;
  d->atom = mol_cif_atom_new(d->n_atoms);
  d->cell = cell_copy_new(c->cell);
  for (i=0; i<c->n_atoms; i++)
    mol_cif_atom_copy(c->atom+i,d->atom+i);
  return(d);
  }

/* copy cif atom data

   a0 - the source atom data structure
   a1 - the new atom data structure */
void mol_cif_atom_copy(struct cif_atom *a0, struct cif_atom *a1) {
  unsigned i;
  str_free(a1->name);
  a1->name = str_copy_new(a0->name);
  a1->num = a0->num;
  for (i=0; i<3; i++)
    a1->coord[i] = a0->coord[i];
  }

/* -------------------------------------------------------------------------- */
