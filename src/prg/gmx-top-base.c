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
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Allocate memory for topology data */
struct gmx_top* gmx_top_new(void) {
  struct gmx_top *d = NULL;
  /* memory allocation */
  d = (struct gmx_top*)malloc(sizeof(struct gmx_top));
  if (!d)
    msg_error("cannot allocate memory for topology data",1);
  /* initialization */
  d->n_atom_types = 0;
  d->n_bonds = 0;
  d->n_angles = 0;
  d->n_dihedrals = 0;
  d->n_impropers = 0;
  d->n_14_pairs = 0;
  d->n_nb_pairs = 0;
  d->n_constraints = 0;
  d->n_settles = 0;
  d->n_igb_parms = 0;
  d->n_cmap = 0;
  d->n_exclusions = 0;
  d->type = NULL;
  d->bond = NULL;
  d->angl = NULL;
  d->dihe = NULL;
  d->impr = NULL;
  d->pair = NULL;
  d->nbpr = NULL;
  d->cnst = NULL;
  d->sttl = NULL;
  d->igbp = NULL;
  d->cmap = NULL;
  d->excl = NULL;
  return(d);
  }

/* Clean memory allocated in topology data structure

   d - the data structure */
void gmx_top_clean(struct gmx_top *d) {
  d->type = gmx_atom_free(d->type,d->n_atom_types);
  d->n_atom_types = 0;
  d->bond = gmx_bond_free(d->bond,d->n_bonds);
  d->n_bonds = 0;
  d->angl = gmx_bond_free(d->angl,d->n_angles);
  d->n_angles = 0;
  d->dihe = gmx_bond_free(d->dihe,d->n_dihedrals);
  d->n_dihedrals = 0;
  d->impr = gmx_bond_free(d->impr,d->n_impropers);
  d->n_impropers = 0;
  d->pair = gmx_bond_free(d->pair,d->n_14_pairs);
  d->n_14_pairs = 0;
  d->nbpr = gmx_bond_free(d->nbpr,d->n_nb_pairs);
  d->n_nb_pairs = 0;
  d->cnst = gmx_bond_free(d->cnst,d->n_constraints);
  d->n_constraints = 0;
  d->sttl = gmx_bond_free(d->sttl,d->n_settles);
  d->n_settles = 0;
  d->igbp = gmx_bond_free(d->igbp,d->n_igb_parms);
  d->n_igb_parms = 0;
  d->cmap = gmx_bond_free(d->cmap,d->n_cmap);
  d->n_cmap = 0;
  d->excl = gmx_bond_free(d->excl,d->n_exclusions);
  d->n_exclusions = 0;
  }

/* Clean memory allocated for topology

   d - the data array */
struct gmx_top* gmx_top_free(struct gmx_top *d) {
  if (d) {
    gmx_top_clean(d);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* Copy topology data from one data structure to another

   d1 - destination data structure
   d2 - source data structure */
void gmx_top_copy(struct gmx_top *d1, struct gmx_top *d2) {
  if (d1 && d2) {
    d1->n_atom_types = d2->n_atom_types;
    d1->n_bonds = d2->n_bonds;
    d1->n_angles = d2->n_angles;
    d1->n_dihedrals = d2->n_dihedrals;
    d1->n_impropers = d2->n_impropers;
    d1->n_14_pairs = d2->n_14_pairs;
    d1->n_nb_pairs = d2->n_nb_pairs;
    d1->n_constraints = d2->n_constraints;
    d1->n_settles = d2->n_settles;
    d1->n_igb_parms = d2->n_igb_parms;
    d1->n_cmap = d2->n_cmap;
    d1->n_exclusions = d2->n_exclusions;
    d1->type = gmx_atom_copy_new(d2->type,d2->n_atom_types);
    d1->bond = gmx_bond_copy_new(d2->bond,d2->n_bonds);
    d1->angl = gmx_bond_copy_new(d2->angl,d2->n_angles);
    d1->dihe = gmx_bond_copy_new(d2->dihe,d2->n_dihedrals);
    d1->impr = gmx_bond_copy_new(d2->impr,d2->n_impropers);
    d1->pair = gmx_bond_copy_new(d2->pair,d2->n_14_pairs);
    d1->nbpr = gmx_bond_copy_new(d2->nbpr,d2->n_nb_pairs);
    d1->cnst = gmx_bond_copy_new(d2->cnst,d2->n_constraints);
    d1->sttl = gmx_bond_copy_new(d2->sttl,d2->n_settles);
    d1->igbp = gmx_bond_copy_new(d2->igbp,d2->n_igb_parms);
    d1->cmap = gmx_bond_copy_new(d2->cmap,d2->n_cmap);
    d1->excl = gmx_bond_copy_new(d2->excl,d2->n_exclusions);
    }
  }

/* Create copy of topology data structure

   d - the data array
   n - length of the array */
struct gmx_top* gmx_top_copy_new(struct gmx_top *d) {
  struct gmx_top *t = NULL;
  if (d) {
    t = gmx_top_new();
    gmx_top_copy(t,d);
    }
  return(t);
  }

/* -------------------------------------------------------------------------- */
