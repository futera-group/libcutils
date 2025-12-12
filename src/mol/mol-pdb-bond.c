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
 
   p - pointer to pdb molecular data struct
   f - residuum ID
   n - number of bonds */
unsigned **mol_pdb_bonds(struct pdb_mol *p, unsigned f, unsigned *n) {
  unsigned i,j,*b,**bn,n1,n2,id=0;
  double d,r;
  struct queue *q;
  q=queue_alloc();
  /* identify bonds */
  for (i=0; i<p->res[f].n_atoms; i++)
    for (j=i+1; j<p->res[f].n_atoms; j++) {
      n1 = atom_num_pdb(p->res[f].atom[i].name);
      n2 = atom_num_pdb(p->res[f].atom[j].name);
      d = dist_r(p->res[f].atom[i].coord,p->res[f].atom[j].coord,3);
      r = 1.25*(atom_radius(n1)+atom_radius(n2));
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
