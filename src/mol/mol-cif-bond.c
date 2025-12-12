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

#include <cmn/matrix.h>
#include <cmn/queue.h>
#include <cmn/vector.h>
#include "mol/atom.h"
#include "mol/distance.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* create bond array with atom indices for the given molecule
 
   c - pointer to xyz molecular data struct
   n - number of bonds */
unsigned **mol_cif_bonds(struct cif_mol *c, unsigned *n) {
 unsigned i,j,k,*b,**bn,id=0;
  double d,r,c1[3],c2[3];
  struct queue *q;
  q = queue_alloc();
  /* identify bonds */
  for (i=0; i<c->n_atoms; i++)
    for (j=i+1; j<c->n_atoms; j++) {
      /* coordinates */
      if (c->cell->transf) {
        mat_fmult_rvec(c1,c->cell->t_c2r_c,c->atom[i].coord,3,3);
        mat_fmult_rvec(c2,c->cell->t_c2r_c,c->atom[j].coord,3,3);
        }
      else {
        vec_fcopy(c1,c->atom[i].coord,3);
        vec_fcopy(c2,c->atom[j].coord,3);
        }
      for (k=0; k<3; k++) {
        c1[k] *= c->cell->side[k];
        c2[k] *= c->cell->side[k];
        }
      /* distance */
      d = dist_r(c1,c2,3);
      r = 1.25*(atom_radius(c->atom[i].num)+atom_radius(c->atom[j].num));
      if (d<r) {
        b = vec_ualloc(2);
        b[0] = i;
        b[1] = j;
        queue_add(q,b);
        }
      }
  /* create bond array */
  (*n) = q->num;
  bn = mat_ualloc(q->num,2);
  while (q->num) {
    b = queue_get(q);
    bn[id][0] = b[0];
    bn[id][1] = b[1];
    vec_ufree(b);
    id++;
    }
  /* clean memory */
  queue_free(q);
  return(bn);
  }

/* -------------------------------------------------------------------------- */
