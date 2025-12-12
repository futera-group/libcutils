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
#include <cmn/stack.h>
#include <cmn/string.h>
#include <cmn/tree.h>
#include <cmn/vector.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert cif molecular data struct to apc struct
 
   c - pointer to cif molecular data struct */
struct apc_mol *mol_cif_apc(struct cif_mol *c) { 
  unsigned i,a,**b,nb,id = 0;
  struct apc_mol *m;
  struct stack *q;
  struct tdata *d;
  struct tree *t;
  /* apc structure */
  m = mol_apc_new();
  b = mol_cif_bonds(c,&nb);
  t = mol_tree_set(b,nb,c->n_atoms,&(m->loop),&(m->n_loops));
  m->title = str_copy_new("Converted from cif file");
  m->file = str_copy_new("prep-mol.res");
  m->resname = str_copy_new("MOL");
  m->n_atoms = c->n_atoms;
  m->atom = mol_apc_atom_new(m->n_atoms);
  /* convert tree to tree marks */
  q = stack_alloc();
  stack_add(q,t->root);
  while (q->num) {
    d = stack_get(q);
    a = ((struct mol_tree_atom*)(d->t_data))->id;
    m->atom[id].name = str_copy_new(c->atom[a].name);
    m->atom[id].type = str_copy_new("x");
    if (c->cell->transf)
      mat_fmult_rvec(m->atom[id].coord,c->cell->t_c2r_c,c->atom[a].coord,3,3);
    else
      vec_fcopy(m->atom[id].coord,c->atom[a].coord,3);
    for (i=0; i<3; i++)
      m->atom[id].coord[i] *= c->cell->side[i];
    m->atom[id].charge = 0.0;
    m->atom[id].tree = ((struct mol_tree_atom*)(d->t_data))->tree;
    m->atom[id].id = a;
    for (i=0; i<d->t_num; i++)
      stack_add(q,d->t_child[i]);
    id++;
    }
  mol_apc_atom_id_set(m);
  /* clean memory */
  mat_ufree(b,nb);
  stack_free(q);
  tree_free(t);
  return(m);
  }

/* -------------------------------------------------------------------------- */
