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

/* allocate new pdb molecule struct */
struct pdb_mol *mol_pdb_new(void) {
  struct pdb_mol *p = NULL;
  /* memory allocation */
  p = (struct pdb_mol*)malloc(sizeof(struct pdb_mol));
  if (!p)
    msg_error("cannot allocate memory for pdb molecular data struct",1);
  /* initialization */
  p->title = NULL;
  p->remark = NULL;
  p->n_res = 0;
  p->res = NULL;
  return(p);
  }

/* allocate memory for array of residui
 
   n - number of residui */
struct pdb_res *mol_pdb_res_new(unsigned n) {
  struct pdb_res *r = NULL;
  unsigned i;
  /* memory allocation */
  r = (struct pdb_res*)malloc(n*sizeof(struct pdb_res));
  if (!r)
    msg_error("cannot allocate memory for array of pdb residui",1);
  /* initilization */
  for (i=0; i<n; i++) {
    r[i].name = NULL;
    r[i].ter = 0;
    r[i].id = 0;
    r[i].n_atoms = 0;
    r[i].atom = NULL;
    }
  return(r);
  }

/* allocate memory for array of atoms
 
   n - number of atoms */
struct pdb_atom *mol_pdb_atom_new(unsigned n) {
  struct pdb_atom *a = NULL;
  unsigned i,j;
  /* memory allocation */
  a = (struct pdb_atom*)malloc(n*sizeof(struct pdb_atom));
  if (!a)
    msg_error("cannot allocate memory for array of pdb atoms",1);
  /* initilization */
  for (i=0; i<n; i++) {
    a[i].name = NULL;
    a[i].id = 0;
    a[i].charge = 0.0;
    for (j=0; j<3; j++)
      a[i].coord[j] = 0.0;
    }
  return(a);
  }

/* allocate memory for new temporary info pdb atom struct */
struct pdb_unit *mol_pdb_unit_new(void) {
  struct pdb_unit *a = NULL;
  unsigned i;
  /* memory allocation */
  a = (struct pdb_unit*)malloc(sizeof(struct pdb_unit));
  if (!a)
    msg_error("cannot allocate memory for new pdb atom",1);
  /* initialization */
  a->ter = 0;
  a->name = NULL;
  a->resname = NULL;
  a->resid = 0;
  a->charge = 0.0;
  for (i=0; i<3; i++)
    a->coord[i] = 0.0;
  return(a);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for pdb molecule struct

   p - pointer to pdb molecular data struct */
void mol_pdb_free(struct pdb_mol *p) {
  if (p) {
    str_free(p->title);
    str_free(p->remark);
    mol_pdb_res_free(p->res,p->n_res);
    free(p);
    }
  }

/* free memory allocated for array of pdb residui 
 
   r - pointer to the array
   n - number of residui */
void mol_pdb_res_free(struct pdb_res *r, unsigned n) {
  unsigned i;
  if (r) {
    for (i=0; i<n; i++) {
      str_free(r[i].name);
      mol_pdb_atom_free(r[i].atom,r[i].n_atoms);
      }
    free(r);
    }
  }

/* free memory allocated for array of pdb atoms 
 
   a - pointer to the array
   n - number of atoms */
void mol_pdb_atom_free(struct pdb_atom *a, unsigned n) {
  unsigned i;
  if (a) {
    for (i=0; i<n; i++)
      str_free(a[i].name);
    free(a);
    }
  }

/* free memory allocated for pbc unit data structure 
 
   a - the pdb unit data structure */
void mol_pdb_unit_free(struct pdb_unit *a) {
  if (a) {
    str_free(a->name);
    str_free(a->resname);
    free(a);
    }
  }

/* -------------------------------------------------------------------------- */
