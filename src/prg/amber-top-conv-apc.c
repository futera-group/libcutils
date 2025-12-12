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
#include <mol/molec.h>
#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

/* convert amber topology and coordinates to apc molecular file format
 
   t - pointer to amber topology struct
   c - pointer to array with coordinates
   r - residuum ID */
struct apc_mol *amber_top_apc(struct amber_top *t, double *c, unsigned r) {
  unsigned i,**b,a,a0,a1,na,nb,n = 0;
  char sym[80];
  short tree;
  struct apc_mol *p;
  struct stack *q;
  struct tdata *d;
  struct tree *s;
  /* first/last atom ID */
  a0 = t->res_pointer[r]-1;
  a1 = (r==(t->pointers[AMBER_POINTER_NRES]-1) ?
        t->pointers[AMBER_POINTER_NATOM]-1 : 
        t->res_pointer[r+1]-1);
  /* number of atoms */
  na = a1-a0+1;
  /* check tree marks */
  tree = 1;
  for (i=a0; i<=a1; i++) {
    str_trim_copy(t->tree_chain_class[i],sym);
    if (str_length(sym)>1) {
      tree = 0;
      break;
      }
    }
  /* conversion to prep */
  p = mol_apc_new();
  p->title = str_sprintf("Residuum #%d: %s",r+1,t->res_names[r]);
  str_lowcase_copy(t->res_names[r],sym);
  str_trim(sym);
  p->file = str_sprintf("prep-%s.res",sym);
  p->resname = str_copy_new(t->res_names[r]);
  p->n_atoms = na;
  p->atom = mol_apc_atom_new(p->n_atoms);
  /* use tree marks from the topology */
  if (tree) {
    for (i=a0; i<=a1; i++) {
      p->atom[i-a0].name = str_copy_new(t->atom_names[i]);
      p->atom[i-a0].type = str_copy_new(t->atom_types[i]);
      p->atom[i-a0].tree = mol_apc_tree_id(t->tree_chain_class[i]);
      vec_fcopy(p->atom[i-a0].coord,c+3*(i-a0),3);
      p->atom[i-a0].charge = t->charge[i]/AMBER_TOP_CHARGE;
      }
    }
  /* create connectivity graph */
  else {
    b = amber_top_bonds_res(t,r,&nb,NULL,AMBER_BOND_ALL);
    s = mol_tree_set(b,nb,na,&(p->loop),&(p->n_loops));
    q = stack_alloc();
    stack_add(q,s->root);
    while (q->num) {
      d = stack_get(q);
      a = ((struct mol_tree_atom*)(d->t_data))->id;
      p->atom[n].name = str_copy_new(t->atom_names[a]);
      p->atom[n].type = str_copy_new(t->atom_types[a]);
      for (i=0; i<3; i++)
        p->atom[n].coord[i] = c[3*a+i];
      p->atom[n].charge = t->charge[a]/AMBER_TOP_CHARGE;
      p->atom[n].tree = ((struct mol_tree_atom*)(d->t_data))->tree;
      p->atom[n].id = a;
      for (i=0; i<d->t_num; i++)
        stack_add(q,d->t_child[i]);
      n++;
      }
    stack_free(q);
    mat_ufree(b,2);
    tree_free(s);
    }
  mol_apc_atom_id_set(p);
  /* clean memory */
  return(p);
  }

/* -------------------------------------------------------------------------- */
