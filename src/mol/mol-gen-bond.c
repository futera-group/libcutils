/******************************************************************************\
 *                                                                            * 
 *  Libcutils - library of C function                                         * 
 *                                                                            *
 *  Version:             3.4                                                  * 
 *  Date:                06/03/2020                                           *
 *                                                                            * 
 *  Author:              Zdenek Futera                                        * 
 *                                                                            * 
 *  Address:             University of South Bohemia                          * 
 *                       Faculty of Science, Institute of Physics             * 
 *                       Branisovska 1760, 370 05 Ceske Budejovice            * 
 *                       Czech Republic                                       * 
 *                                                                            * 
 *  E-Mail:              zfutera@prf.jcu.cz                                   * 
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
 
   x - molecular data
   n - number of bonds (output) */
unsigned **mol_gen_bonds(struct gen_mol *x, unsigned *n) {
  unsigned i,j,*b,**bn,id = 0;
  double d,r;
  struct queue *q;
  q = queue_alloc();
  /* identify bonds */
  for (i=0; i<x->n_atoms; i++)
    for (j=i+1; j<x->n_atoms; j++) {
      d = dist_r(x->atom[i].coord,x->atom[j].coord,3);
      r = 1.25*(atom_radius(x->atom[i].num)+atom_radius(x->atom[j].num));
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
