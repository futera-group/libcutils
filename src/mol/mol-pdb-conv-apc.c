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
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert pdb molecular data struct to apc struct
 
   p  - pointer to pdb molecular data struct
   id - residuum ID */
struct apc_mol *mol_pdb_apc(struct pdb_mol *p, unsigned id) {
  char s0[1024],s1[1024];
  unsigned i,**b,a,nb,n = 0;
  struct apc_mol *r;
  struct stack *q;
  struct tdata *d;
  struct tree *t;
  /* apc structure */
  r = mol_apc_new();
  b = mol_pdb_bonds(p,id,&nb);
  t = mol_tree_set(b,nb,p->res[id].n_atoms,&(r->loop),&(r->n_loops));
  sprintf(s0,"Residuum #%d: %s",id+1,p->res[id].name);
  r->title = str_copy_new(s0);
  str_lowcase_copy(p->res[id].name,s0);
  sprintf(s1,"prep-%s.res",s0);
  r->file = str_copy_new(s1);
  r->resname = str_copy_new(p->res[id].name);
  r->n_atoms = p->res[id].n_atoms;
  r->atom = mol_apc_atom_new(r->n_atoms);
  /* convert tree to tree marks */
  q = stack_alloc();
  stack_add(q,t->root);
  while (q->num) {
    d = stack_get(q);
    a = ((struct mol_tree_atom*)(d->t_data))->id;
    r->atom[n].name = str_copy_new(p->res[id].atom[a].name);
    r->atom[n].type = str_copy_new("x");
    for (i=0; i<3; i++)
      r->atom[n].coord[i] = p->res[id].atom[a].coord[i];
    r->atom[n].charge = p->res[id].atom[a].charge;
    r->atom[n].tree = ((struct mol_tree_atom*)(d->t_data))->tree;
    r->atom[n].id = a;
    for (i=0; i<d->t_num; i++)
      stack_add(q,d->t_child[i]);
    n++;
    }
  mol_apc_atom_id_set(r);
  /* clean memory */
  stack_free(q);
  mat_ufree(b,nb);
  tree_free(t);
  return(r);
  }

/* -------------------------------------------------------------------------- */
