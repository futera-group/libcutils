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

/* convert xyz molecular data struct to apc struct
 
   x - pointer to xyz molecular data struct */
struct apc_mol *mol_xyz_apc(struct xyz_mol *x) { 
  unsigned i,a,id,**b,nb,at[120];
  char sym[128];
  struct apc_mol *m;
  struct stack *q;
  struct tdata *d;
  struct tree *t;
  /* initialization */
  vec_uset(at,0,120);
  id = 0;
  /* apc structure */
  m = mol_apc_new();
  b = mol_xyz_bonds(x,&nb);
  t = mol_tree_set(b,nb,x->n_atoms,&(m->loop),&(m->n_loops));
  m->title = str_copy_new("Converted from xyz file");
  m->file = str_copy_new("prep-mol.res");
  m->resname = str_copy_new("MOL");
  m->n_atoms = x->n_atoms;
  m->atom = mol_apc_atom_new(m->n_atoms);
  /* convert tree to tree marks */
  q = stack_alloc();
  stack_add(q,t->root);
  while (q->num) {
    d = stack_get(q);
    a = ((struct mol_tree_atom*)(d->t_data))->id;
    if (x->atom[a].name)
      m->atom[id].name = str_copy_new(x->atom[a].name);
    else {
      mol_pdb_atom_name_id(sym,x->atom[a].num,at[x->atom[a].num]++);
      m->atom[id].name = str_copy_new(sym);
      }
    m->atom[id].type = str_copy_new("x");
    for (i=0; i<3; i++)
      m->atom[id].coord[i] = x->atom[a].coord[i];
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
