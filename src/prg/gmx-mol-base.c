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
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Allocate memory for array of molecules

   n - number of molecules */
struct gmx_mol* gmx_mol_new(unsigned n) {
  unsigned i;
  struct gmx_mol *d = NULL;
  /* memory allocation */
  d = (struct gmx_mol*)malloc(n*sizeof(struct gmx_mol));
  if (!d)
    msg_error("cannot allocate memory for molecular data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].file = NULL;
    d[i].name = NULL;
    d[i].gen_itp = 0;
    d[i].n_excl = 0;
    d[i].n_atoms = 0;
    d[i].n_resids = 0;
    d[i].res = NULL;
    d[i].top = gmx_top_new();
    }
  return(d);
  }

/* Clean memory allocated in molecule structure

   d - the data structure */
void gmx_mol_clean(struct gmx_mol *d) {
  d->file = str_free(d->file);
  d->name = str_free(d->name);
  d->n_excl = 0;
  d->n_atoms = 0;
  d->res = gmx_res_free(d->res,d->n_resids);
  d->n_resids = 0;
  d->top = gmx_top_free(d->top);
  }

/* Clean memory allocated for molecule array

   d - the data array
   n - length of the array */
struct gmx_mol* gmx_mol_free(struct gmx_mol *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++)
      gmx_mol_clean(d+i);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* Copy molecule data from one data structure to another

   d1 - destination data structure
   d2 - source data structure */
void gmx_mol_copy(struct gmx_mol *d1, struct gmx_mol *d2) {
  if (d1 && d2) {
    d1->file = str_copy_new(d2->file);
    d1->name = str_copy_new(d2->name);
    d1->gen_itp = d2->gen_itp;
    d1->n_excl = d2->n_excl;
    d1->n_atoms = d2->n_atoms;
    d1->n_resids = d2->n_resids;
    d1->res = gmx_res_copy_new(d2->res,d2->n_resids);
    d1->top = gmx_top_copy_new(d2->top);
    }
  }

/* Create copy of molecule data structure array

   d - the data array
   n - length of the array */
struct gmx_mol* gmx_mol_copy_new(struct gmx_mol *d, unsigned n) {
  unsigned i;
  struct gmx_mol *t = NULL;
  if (d) {
    t = gmx_mol_new(n);
    for (i=0; i<n; i++)
      gmx_mol_copy(t+i,d+i);
    }
  return(t);
  }

/* -------------------------------------------------------------------------- */

/* List searching function comparing molecule file names

   d1,d2 - molecular data structures */
short gmx_mol_cmp_file_wrp_search(void *d1, void *d2) {
  if (!d1 || !d2)
    return(0);
  return(str_compare(
    ((struct gmx_mol*)d1)->file,
    ((struct gmx_mol*)d2)->file));
  }

/* List sorting function comparing molecule file names

   d1,d2 - molecular data structures */
short gmx_mol_cmp_file_wrp_sort(void *d1, void *d2) {
  if (!d1 || !d2)
    return(0);
  return(strcmp(
    ((struct gmx_mol*)d1)->file,
    ((struct gmx_mol*)d2)->file));
  }

/* Compare two molecular data structures
 
   m1,m2 - the data structures */
short gmx_mol_compare(struct gmx_mol *m1, struct gmx_mol *m2) {
  unsigned i;
  /* name */
  if (!str_compare(m1->name,m2->name))
    return(0);
  /* number of atoms */
  if (m1->n_atoms!=m2->n_atoms)
    return(0);
  /* number of residues */
  if (m1->n_resids!=m2->n_resids)
    return(0);
  /* residues */
  for (i=0; i<m1->n_resids; i++)
    if (!gmx_res_compare(m1->res+i,m2->res+i))
      return(0); 
  return(1);
  }

/* -------------------------------------------------------------------------- */
