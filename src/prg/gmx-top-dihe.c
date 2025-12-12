
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
#include <cmn/queue.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Sort save dihedral types / definition to proper / improper

   d - gromacs topology data structure */
void gmx_top_dihe_prop_impr(struct gmx_top *d) {
  unsigned i,id,np,nm;
  struct gmx_bond *bp,*bm;
  struct queue *q1p,*q1m,*q2p,*q2m;
  /* allocate queues for data IDs */
  q1p = queue_alloc();
  q1m = queue_alloc();
  q2p = queue_alloc();
  q2m = queue_alloc();
  /* collect data IDs */
  for (i=0; i<d->n_dihedrals; i++)
    if (d->dihe[i].type==2)
      queue_uadd(q1m,i);
    else 
      queue_uadd(q1p,i);
  for (i=0; i<d->n_impropers; i++)
    if (d->impr[i].type==2)
      queue_uadd(q2m,i);
    else 
      queue_uadd(q2p,i);
  /* create new data arrays */
  np = q1p->num + q2p->num;
  bp = gmx_bond_new(np);
  id = 0;
  while (q1p->num) {
    queue_uget(q1p,&i);
    gmx_bond_copy(bp+id,d->dihe+i);
    id++;
    }
  while (q2p->num) {
    queue_uget(q2p,&i);
    gmx_bond_copy(bp+id,d->impr+i);
    id++;
    }
  nm = q1m->num + q2m->num;
  bm = gmx_bond_new(nm);
  id = 0;
  while (q1m->num) {
    queue_uget(q1m,&i);
    gmx_bond_copy(bm+id,d->dihe+i);
    id++;
    }
  while (q2m->num) {
    queue_uget(q2m,&i);
    gmx_bond_copy(bm+id,d->impr+i);
    id++;
    }
  /* save sorted data */
  d->dihe = gmx_bond_free(d->dihe,d->n_dihedrals);
  d->n_dihedrals = np;
  d->dihe = bp;
  d->impr = gmx_bond_free(d->impr,d->n_impropers);
  d->n_impropers = nm;
  d->impr = bm;
  /* clean memory */
  queue_free(q1p);
  queue_free(q1m);
  queue_free(q2p);
  queue_free(q2m);
  }

/* -------------------------------------------------------------------------- */
