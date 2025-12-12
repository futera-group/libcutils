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

#include <stdio.h>
#include <cmn/matrix.h>
#include <cmn/queue.h>
#include <cmn/string.h>
#include <cmn/vector.h>
#include "mol/atom.h"
#include "mol/distance.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* set empirical formula of molecular system in acf format
 
   m - pointer to acf molecular data struct */
void mol_acf_set_formula(struct acf_mol *m) {
  char *p,s[256] = "\0";
  unsigned i,at[120];
  short first;
  /* count atom types */
  vec_uset(at,0,120);
  for (i=0; i<m->n_atoms; i++)
    at[m->atom[i].num]++;
  /* create formula */
  for (i=0, first=1; i<120; i++)
    if (at[i]) {
      if (first) {
        sprintf(s,"%s%d",atom_name(i),at[i]);
        first = 0;
        }
      else
        sprintf(s," %s%d",atom_name(i),at[i]);
      p = str_merge_new(m->form,s);
      str_free(m->form);
      m->form = p;
      }
  }

/* calculate bonds for molecule in acf format
 
   c - pointer to acf molecular data struct */
void mol_acf_set_bonds(struct acf_mol *c) {
  unsigned i,j,*b,id = 0;
  double d,r;
  struct queue *q;
  /* initialization */
  c->bond = mat_ufree(c->bond,c->n_bonds);
  c->n_bonds=0;
  q = queue_alloc();
  /* identify bonds */
  for (i=0; i<c->n_atoms; i++)
    for (j=i+1; j<c->n_atoms; j++) {
      d = dist_r(c->atom[i].coord,c->atom[j].coord,3);
      r = 1.25*(atom_radius(c->atom[i].num)+atom_radius(c->atom[j].num));
      if (d<r) {
        b = vec_ualloc(2);
        b[0] = i;
        b[1] = j;
        queue_add(q,b);
        }
      }
  /* create bond array */
  c->n_bonds = q->num;
  c->bond = mat_ualloc(c->n_bonds,2);
  while (q->num) {
    b = queue_get(q);
    c->bond[id][0] = b[0];
    c->bond[id][1] = b[1];
    vec_ufree(b);
    id++;
    }
  /* clean memory */
  queue_free(q);
  }

/* -------------------------------------------------------------------------- */
