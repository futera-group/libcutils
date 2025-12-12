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
#include <cmn/vector.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Allocate memory for gromacs data structure */
struct gmx_dat* gmx_dat_new(void) {
  struct gmx_dat *d = NULL;
  /* memory allocation */
  d = (struct gmx_dat*)malloc(sizeof(struct gmx_dat));
  if (!d)
    msg_error("cannot allocate memory for gromacs data",1);
  /* initialization */
  d->title = NULL;
  d->n_atoms = 0;
  d->n_resids = 0;
  d->n_mols = 0;
  d->n_frags = 0;
  d->box = NULL;
  d->crd = NULL;
  d->ff = gmx_ff_new();
  d->mol = NULL;
  d->frag = NULL;
  return(d);
  }

/* Clean memory allocated in gromacs data structure

   d - the data structure */
void gmx_dat_clean(struct gmx_dat *d) {
  d->title = str_free(d->title);
  d->n_resids = 0; 
  d->box = vec_ffree(d->box);
  d->crd = mat_ffree(d->crd,d->n_atoms);
  d->n_atoms = 0;
  d->ff = gmx_ff_free(d->ff);
  d->mol = gmx_mol_free(d->mol,d->n_mols);
  d->n_mols = 0;
  d->frag = gmx_frag_free(d->frag,d->n_frags);
  d->n_frags = 0;
  }

/* Clean memory allocated for gromacs data structure

   d - the data array */
struct gmx_dat* gmx_dat_free(struct gmx_dat *d) {
  if (d) {
    gmx_dat_clean(d);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* Copy gromacs data from one data structure to another

   d1 - destination data structure
   d2 - source data structure */
void gmx_dat_copy(struct gmx_dat *d1, struct gmx_dat *d2) {
  if (d1 && d2) {
    d1->title = str_copy_new(d2->title);
    d1->n_atoms = d2->n_atoms;
    d1->n_resids = d2->n_resids;
    d1->n_mols = d2->n_mols;
    d1->n_frags = d2->n_frags;
    d1->box = vec_fcopy_new(d2->box,3);
    d1->crd = mat_fcopy_new(d2->crd,d2->n_atoms,3);
    d1->ff = gmx_ff_copy_new(d2->ff);
    d1->mol = gmx_mol_copy_new(d2->mol,d2->n_mols);
    d1->frag = gmx_frag_copy_new(d2->frag,d2->n_frags);
    }
  }

/* Create copy of gromacs data structure

   d - the data array
   n - length of the array */
struct gmx_dat* gmx_dat_copy_new(struct gmx_dat *d) {
  struct gmx_dat *t = NULL;
  if (d) {
    t = gmx_dat_new();
    gmx_dat_copy(t,d);
    }
  return(t);
  }

/* -------------------------------------------------------------------------- */

/* Update number of atoms and residues in the system by fragment definitions
 
   d - gromacs data structure */
void gmx_dat_update_nums(struct gmx_dat *d) {
  unsigned i;
  d->n_atoms = 0;
  d->n_resids = 0;
  for (i=0; i<d->n_frags; i++) {
    d->n_resids += d->frag[i].n_rep*d->frag[i].mol->n_resids;
    d->n_atoms += d->frag[i].n_rep*d->frag[i].mol->n_atoms;
    }
  }

/* -------------------------------------------------------------------------- */

/* Return system ID of specified atom in the system
 
   g  - gromacs data structure
   ig - fragment ID
   ip - fragment duplicate
   ir - residuum ID
   ia - atom ID */
unsigned gmx_dat_atom_get_id(struct gmx_dat *d, unsigned ig, unsigned ip,
  unsigned ir, unsigned ia) {
  unsigned i,n;
  for (i=0,n=0; i<ig; i++)
    n += d->frag[i].n_rep*d->frag[i].mol->n_atoms;
  n += ip*d->frag[i].mol->n_atoms;
  for (i=0; i<ir; i++)
    n += d->frag[ig].mol->res[i].n_atoms;
  n += ia;
  return(n);
  }

/* -------------------------------------------------------------------------- */
