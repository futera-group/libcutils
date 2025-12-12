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

/* allocate new zmt molecule struct */
struct zmt_mol *mol_zmt_new(void) {
  struct zmt_mol *z = NULL;
  /* memory allocation */
  z = (struct zmt_mol*)malloc(sizeof(struct zmt_mol));
  if (!z)
    msg_error("cannot allocate memory for zmt molecular data struct",1);
  /* initialization */
  z->n_atoms = 0;
  z->atom = NULL;
  return(z);
  }

/* allocate new zmt atomic struct
 
   n - number of atoms in the array */
struct zmt_atom *mol_zmt_atom_new(unsigned n) {
  unsigned i;
  struct zmt_atom *z = NULL;
  /* memory allocation */ 
  z = (struct zmt_atom*)malloc(n*sizeof(struct zmt_atom));
  if (!z)
    msg_error("cannot allocate memory for zmt atomic data struct",1);
  /* initialization */
  for (i=0; i<n; i++) {
    z[i].name = NULL;
    z[i].num = 0;
    z[i].bond_id = 0;
    z[i].angle_id = 0;
    z[i].dihed_id = 0;
    z[i].bond_val = 0.0;
    z[i].angle_val = 0.0;
    z[i].dihed_val = 0.0;
    }
  return(z);
  }

/* allocate new zmt line data structure */
struct zmt_line *mol_zmt_line_new(void) {
  struct zmt_line *z = NULL;
  /* memory allocation */
  z = (struct zmt_line*)malloc(sizeof(struct zmt_line));
  if (!z)
    msg_error("cannot allocate memory for zmt internal atomic data struct",1);
  /* initialization */
  z->num = 0;
  z->bond_id = 0;
  z->bond_sym = NULL;
  z->angle_id = 0;
  z->angle_sym = NULL;
  z->dihed_id = 0;
  z->dihed_sym = NULL;
  return(z);
  }

/* allocate memory for a new zmt internal-coordinate data structure */
struct zmt_icrd *mol_zmt_icrd_new(void) {
  /* memory allocation */
  struct zmt_icrd *z = NULL;
  z = (struct zmt_icrd*)malloc(sizeof(struct zmt_icrd));
  if (!z)
    msg_error("cannot allocate memory for zmt internal-coordinate data",1);
  /* initialization */
  z->sym = NULL;
  z->val = 0.0;
  return(z);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for zmt molecule struct

   z - pointer to zmt molecular data struct */
void mol_zmt_free(struct zmt_mol *z) {
  if (z) {
    mol_zmt_atom_free(z->atom);
    free(z);
    }
  }

/* free memory allocated for array of zmt atom data structures

   a - array of zmt atom data structures */
void mol_zmt_atom_free(struct zmt_atom *a) {
  if (a) {
    a->name = str_free(a->name);
    free(a);
    }
  }

/* free memory allocated for zmt line data structure

   z - the zmt line data structure */
void mol_zmt_line_free(struct zmt_line *z) {
  if (z) {
    str_free(z->bond_sym);
    str_free(z->angle_sym);
    str_free(z->dihed_sym);
    free(z);
    }
  }

/* free memory allocated for zmt internal-coordinae data structure

   z - the zmt internal-coordinate data structure */
void mol_zmt_icrd_free(struct zmt_icrd *z) {
  if (z) {
    str_free(z->sym);
    free(z);
    }
  }

/* -------------------------------------------------------------------------- */

/* copy zmt atom data

   a0 - the source atom data structure
   a1 - the new atom data structure */
void mol_zmt_atom_copy(struct zmt_atom *a0, struct zmt_atom *a1) {
  a1->name = str_free(a1->name);
  a1->name = str_copy_new(a0->name);
  a1->num = a0->num;
  a1->bond_id = a0->bond_id;
  a1->bond_val = a0->bond_val;
  a1->angle_id = a0->angle_id;
  a1->angle_val = a0->angle_val;
  a1->dihed_id = a0->dihed_id;
  a1->dihed_val = a0->dihed_val;
  }

/* -------------------------------------------------------------------------- */
