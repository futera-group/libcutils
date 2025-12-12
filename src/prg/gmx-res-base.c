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

/* Allocate memory for array of residues

   n - number of residues */
struct gmx_res* gmx_res_new(unsigned n) {
  unsigned i;
  struct gmx_res *d = NULL;
  /* memory allocation */
  d = (struct gmx_res*)malloc(n*sizeof(struct gmx_res));
  if (!d)
    msg_error("cannot allocate memory for residuum data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].name = NULL;
    d[i].id = 0;
    d[i].n_atoms = 0;
    d[i].atom = NULL;
    d[i].top = gmx_top_new();
    }
  return(d);
  }

/* Clean memory allocated in residuum structure

   d - the data structure */
void gmx_res_clean(struct gmx_res *d) {
  d->name = str_free(d->name);
  d->id = 0;
  d->atom = gmx_atom_free(d->atom,d->n_atoms);
  d->n_atoms = 0;
  d->top = gmx_top_free(d->top);
  }

/* Clean memory allocated for residuum array

   d - the data array
   n - length of the array */
struct gmx_res* gmx_res_free(struct gmx_res *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++)
      gmx_res_clean(d+i);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* Copy residuum data from one data structure to another

   d1 - destination data structure
   d2 - source data structure */
void gmx_res_copy(struct gmx_res *d1, struct gmx_res *d2) {
  if (d1 && d2) {
    d1->name = str_copy_new(d2->name);
    d1->id = d2->id;
    d1->n_atoms = d2->n_atoms;
    d1->atom = gmx_atom_copy_new(d2->atom,d2->n_atoms);
    d1->top = gmx_top_copy_new(d2->top);
    }
  }

/* Create copy of residuum data structure array

   d - the data array
   n - length of the array */
struct gmx_res* gmx_res_copy_new(struct gmx_res *d, unsigned n) {
  unsigned i;
  struct gmx_res *t = NULL;
  if (d) {
    t = gmx_res_new(n);
    for (i=0; i<n; i++)
      gmx_res_copy(t+i,d+i);
    }
  return(t);
  }

/* -------------------------------------------------------------------------- */

/* Compare two residuum data structures
 
   r1,r2 - the data structures */
short gmx_res_compare(struct gmx_res *r1, struct gmx_res *r2) {
  unsigned i;
  /* name */
  if (!str_compare(r1->name,r2->name))
    return(0);
  /* number of atoms */
  if (r1->n_atoms!=r2->n_atoms)
    return(0);
  /* atoms */
  for (i=0; i<r1->n_atoms; i++)
    if (!gmx_atom_compare(r1->atom+i,r2->atom+i))
      return(0);
  return(1);
  }

/* -------------------------------------------------------------------------- */
