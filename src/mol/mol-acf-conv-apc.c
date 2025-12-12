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

#include <cmn/stack.h>
#include <cmn/string.h>
#include <cmn/tree.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert acf molecular data struct to apc struct
 
   p - pointer to amber ac data struct */
struct apc_mol *mol_acf_apc(struct acf_mol *p) { 
  char s0[1024],s1[1024];
  unsigned i,a,id = 0;
  struct apc_mol *m;
  struct stack *q;
  struct tdata *d;
  struct tree *t;
  /* apc structure */
  m = mol_apc_new();
  mol_acf_connect(p);
  t = mol_tree_set(p->bond,p->n_bonds,p->n_atoms,&(m->loop),&(m->n_loops));
  sprintf(s0,"%s: %s",p->name,p->form);
  m->title = str_copy_new(s0);
  str_lowcase_copy(p->name,s0);
  sprintf(s1,"prep-%s.res",s0);
  m->file = str_copy_new(s1);
  m->resname = str_copy_new(p->name);
  m->n_atoms = p->n_atoms;
  m->atom = mol_apc_atom_new(m->n_atoms);
  /* convert tree to tree marks */
  q = stack_alloc();
  stack_add(q,t->root);
  while (q->num) {
    d = stack_get(q);
    a = ((struct mol_tree_atom*)(d->t_data))->id;
    m->atom[id].name = str_copy_new(p->atom[a].name);
    m->atom[id].type = str_copy_new(p->atom[a].type);
    for (i=0; i<3; i++)
      m->atom[id].coord[i] = p->atom[a].coord[i];
    m->atom[id].charge = p->atom[a].charge;
    m->atom[id].tree = ((struct mol_tree_atom*)(d->t_data))->tree;
    m->atom[id].id = a;
    for (i=0; i<d->t_num; i++)
      stack_add(q,d->t_child[i]);
    id++;
    }
  mol_apc_atom_id_set(m);
  /* clean memory */
  stack_free(q);
  tree_free(t);
  return(m);
  }

/* -------------------------------------------------------------------------- */
