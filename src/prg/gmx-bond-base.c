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
#include <cmn/vector.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Allocate memory for array of bonds

   n - number of bonds */
struct gmx_bond* gmx_bond_new(unsigned n) {
  unsigned i;
  struct gmx_bond *d = NULL;
  /* memory allocation */
  d = (struct gmx_bond*)malloc(n*sizeof(struct gmx_bond));
  if (!d)
    msg_error("cannot allocate memory for bonding data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].type = 0;
    d[i].n_atoms = 0;
    d[i].n_parms = 0;
    d[i].name = NULL;
    d[i].ia = NULL;
    d[i].ir = NULL;
    d[i].parm = NULL;
    }
  return(d);
  }

/* Clean memory allocated in bond data structure

   d - the data structure */
void gmx_bond_clean(struct gmx_bond *d) {
  d->type = 0;
  d->n_parms = 0;
  d->name = vec_sfree(d->name,d->n_atoms);
  d->n_atoms = 0;
  vec_ufree(d->ia);
  vec_ufree(d->ir);
  vec_ffree(d->parm);
  }

/* Clean memory allocated for bond array

   d - the data array
   n - length of the array */
struct gmx_bond* gmx_bond_free(struct gmx_bond *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++)
      gmx_bond_clean(d+i);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* Copy data from one bond data structure to another

   d1 - destination data structure
   d2 - source data structure */
void gmx_bond_copy(struct gmx_bond *d1, struct gmx_bond *d2) {
  unsigned i;
  if (d1 && d2) {
    d1->type = d2->type;
    d1->n_atoms = d2->n_atoms;
    d1->n_parms = d2->n_parms;
    if (d2->name) {
      d1->name = vec_talloc(sizeof(char*),d2->n_atoms);
      for (i=0; i<d2->n_atoms; i++)
        d1->name[i] = str_copy_new(d2->name[i]);
      }
    d1->ia = vec_ucopy_new(d2->ia,d2->n_atoms);
    d1->ir = vec_ucopy_new(d2->ir,d2->n_atoms);
    d1->parm = vec_fcopy_new(d2->parm,d2->n_parms);
    }
  }

/* Create copy of bonding data structure array

   d - the data array
   n - length of the array */
struct gmx_bond* gmx_bond_copy_new(struct gmx_bond *d, unsigned n) {
  unsigned i;
  struct gmx_bond *t = NULL;
  if (d) {
    t = gmx_bond_new(n);
    for (i=0; i<n; i++)
      gmx_bond_copy(t+i,d+i);
    }
  return(t);
  }

/* -------------------------------------------------------------------------- */

/* Add new terms to bonding data array 

   b0 - the bonding data array
   n0 - number of terms
   b1 - new bonding terms
   n1 - number of new terms */
void gmx_bond_add(struct gmx_bond **b0, unsigned *n0,
  struct gmx_bond *b1, unsigned n1) {
  unsigned i,nt;
  struct gmx_bond *t;
  /* extended bonding array */
  nt = (*n0) + n1;
  t = gmx_bond_new(nt);
  for (i=0; i<(*n0); i++)
    gmx_bond_copy(t+i,(*b0)+i);
  for (i=0; i<n1; i++)
    gmx_bond_copy(t+(*n0)+i,b1+i);
  /* replate the array */
  (*b0) = gmx_bond_free(*b0,*n0);
  (*n0) = nt;
  (*b0) = t;
  }

/* Delete specified term from bonding data array 

   b  - the bonding data array
   id - ID of the term to delete
   n  - number of terms */
void gmx_bond_del_one(struct gmx_bond **b, unsigned id, unsigned *n) {
  unsigned i;
  struct gmx_bond *t;
  /* sanity check */
  if (!(*b))
    msg_error("attempt to delete term from empty bonding array",1);
  if (id>=(*n))
    msg_error_f("ID of bonding term to delete out of range"
     " (%d/%d)",1,id+1,*n);
  /* remove the term */
  t = gmx_bond_new((*n)-1);
  for (i=0; i<id; i++)
    gmx_bond_copy(t+i,(*b)+i);
  for (i=id+1; i<(*n); i++)
    gmx_bond_copy(t+i-1,(*b)+i);
  /* replate the array */
  (*b) = gmx_bond_free(*b,*n);
  (*n) = (*n)-1;
  (*b) = t;
  }

/* Delete specified terms from bonding data array 

   b  - the bonding data array
   id - ID of the term to delete
   nd - number of the IDs 
   n  - number of terms */
void gmx_bond_del(struct gmx_bond **b, unsigned *id, unsigned nd, unsigned *n) {
  unsigned i,it,nt;
  short *a;
  struct gmx_bond *t;
  if (!nd)
    return;
  /* sanity check */
  if (!(*b))
    msg_error("attempt to delete term from empty bonding array",1);
  for (i=0; i<nd; i++)
    if (id[i]>=(*n))
      msg_error_f("ID of bonding term to delete out of range"
        " (%d/%d)",1,id[i]+1,*n);
  /* array of delete/stay indicators */
  a = vec_sialloc(*n);
  vec_siset(a,1,*n);
  for (i=0; i<nd; i++)
    a[id[i]] = 0;
  /* number of terms after deletion */
  for (i=0,nt=0; i<(*n); i++)
    if (a[i])
      nt++;
  /* delete all */
  if (nt==0) {
    (*b) = gmx_bond_free(*b,*n);
    (*n) = 0;
    }
  /* delete selected */
  else if (nt<(*n)) {
    t = gmx_bond_new(nt);
    for (i=0,it=0; i<(*n); i++)
      if (a[i])
        gmx_bond_copy(t+it++,(*b)+i);
    (*b) = gmx_bond_free(*b,*n);
    (*n) = nt;
    (*b) = t;
    }
  /* clean memory */
  vec_sifree(a);
  }

/* -------------------------------------------------------------------------- */
