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

#include <cmn/message.h>
#include <cmn/queue.h>
#include <cmn/vector.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Return pointer to specified bonding term 

   d  - topology data structure
   id - internal ID of the bonding term 
   n  - number of the terms (output,optional) */
struct gmx_bond** gmx_top_bond_get(struct gmx_top *d, short id, unsigned **n) {
  switch (id) {
    /* bonds */
    case GMX_BTERM_BOND:
      if (n) (*n) = &(d->n_bonds);
      return(&(d->bond));
    /* angles */
    case GMX_BTERM_ANGL:
      if (n) (*n) = &(d->n_angles);
      return(&(d->angl));
    /* dihedrals */
    case GMX_BTERM_DIHE:
      if (n) (*n) = &(d->n_dihedrals);
      return(&(d->dihe));
    /* impropers */
    case GMX_BTERM_IMPR:
      if (n) (*n) = &(d->n_impropers);
      return(&(d->impr));
    /* 1-4 pairs */
    case GMX_BTERM_PAIR:
      if (n) (*n) = &(d->n_14_pairs);
      return(&(d->pair));
    /* other non-bonded pairs */
    case GMX_BTERM_NBPR:
      if (n) (*n) = &(d->n_nb_pairs);
      return(&(d->nbpr));
    /* constraints */
    case GMX_BTERM_CNST:
      if (n) (*n) = &(d->n_constraints);
      return(&(d->cnst));
    /* settles */
    case GMX_BTERM_STTL:
      if (n) (*n) = &(d->n_settles);
      return(&(d->sttl));
    /* cmap corrections */
    case GMX_BTERM_CMAP:
      if (n) (*n) = &(d->n_cmap);
      return(&(d->cmap));
    /* exclusions */
    case GMX_BTERM_EXCL:
      if (n) (*n) = &(d->n_exclusions);
      return(&(d->excl));
    }
  msg_error_f("unknown bonding-term ID %d",1,id);
  return(NULL);
  }

/* -------------------------------------------------------------------------- */

/* Remove all bonding terms involving specified atom
 
   d   - topology data structure
   ir  - residuum ID 
   ia  - ID of the atom
   upd - update atomic pointers */
void gmx_top_bond_del(struct gmx_top *d, unsigned ir, unsigned ia, short upd) {
  short term[10] = {GMX_BTERM_BOND,GMX_BTERM_ANGL,GMX_BTERM_DIHE,
    GMX_BTERM_IMPR,GMX_BTERM_PAIR,GMX_BTERM_NBPR,GMX_BTERM_CNST,
    GMX_BTERM_STTL,GMX_BTERM_CMAP,GMX_BTERM_EXCL};
  unsigned i,j,k,*t,nt,*n,n_terms = 10;
  struct gmx_bond **b;
  struct queue *q;
  /* storage for bond IDs */
  q = queue_alloc();
  /* collect the IDs */
  for (i=0; i<n_terms; i++) {
    b = gmx_top_bond_get(d,term[i],&n);
    for (j=0; j<(*n); j++) {
      for (k=0; k<(*b)[j].n_atoms; k++) {
        if ((*b)[j].ir[k]==ir && (*b)[j].ia[k]==ia) {
          queue_uadd(q,j);
          break;
          }
        }
      }
    /* delete selected terms */
    nt = 0;
    t = vec_ualloc(q->num);
    while (q->num)
      queue_uget(q,&(t[nt++]));
    gmx_bond_del(b,t,nt,n);
    t = vec_ufree(t);
    /* update atomic pointers */
    if (upd) {
      for (j=0; j<(*n); j++)
        for (k=0; k<(*b)[j].n_atoms; k++) {
          if ((*b)[j].ir[k]==ir && (*b)[j].ia[k]>ia)
            (*b)[j].ia[k]--;
          }
      }
    }
  /* clean memory */
  queue_free(q);
  }   

/* -------------------------------------------------------------------------- */
