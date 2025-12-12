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
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* allocate new acf molecule data */
struct acf_mol *mol_acf_new(void) {
  struct acf_mol *c = NULL;
  /* memory allocation */
  c = (struct acf_mol*)malloc(sizeof(struct acf_mol));
  if (!c)
    msg_error("cannot allocate memory for acf molecular data struct",1);
  /* initialization */
  c->form = NULL;
  c->name = NULL;
  c->charge = 0.0;
  c->n_atoms = 0;
  c->n_bonds = 0;
  c->bond = NULL;
  c->atom = NULL;
  return(c);
  }

/* allocate new acf atom data

   n - number of atoms */
struct acf_atom *mol_acf_atom_new(unsigned n) {
  unsigned i,j;
  struct acf_atom *c = NULL;
  /* memory allocation */
  c = (struct acf_atom*)malloc(n*sizeof(struct acf_atom));
  if (!c)
    msg_error("cannot allocate memory for acf atom data struct",1);
  /* initialization */
  for (i=0; i<n; i++) {
    c[i].name = NULL;
    c[i].type = NULL;
    c[i].num = 0;
    c[i].res = 0;
    c[i].charge = 0.0;
    for (j=0; j<3; j++)
      c[i].coord[j] = 0.0;
    }
  return(c);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for acf molecule data
 
   c - the acf molecular data */
void mol_acf_free(struct acf_mol *c) {
  if (c) {
    c->bond = mat_ufree(c->bond,c->n_bonds);
    mol_acf_atom_free(c->atom,c->n_atoms);
    free(c);
    }
  }

/* free memory allocated for acf atom data
 
   c - the acf atom data
   n - number of atoms */
void mol_acf_atom_free(struct acf_atom *c, unsigned n) {
  unsigned i;
  if (c) {
    for (i=0; i<n; i++) {
      str_free(c[i].name);
      str_free(c[i].type);
      }
    free(c);
    }
  }

/* -------------------------------------------------------------------------- */

/* copy acf atom data

   a0 - the source atom data structure
   a1 - the new atom data structure */
void mol_acf_atom_copy(struct acf_atom *a0, struct acf_atom *a1) {
  unsigned i;
  str_free(a1->name);
  a1->name = str_copy_new(a0->name);
  str_free(a1->type);
  a1->type = str_copy_new(a0->type);
  a1->num = a0->num;
  a1->res = a0->res;
  a1->charge = a0->charge;
  for (i=0; i<3; i++)
    a1->coord[i] = a0->coord[i];
  }

/* -------------------------------------------------------------------------- */
