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

#include <stdlib.h>
#include <cmn/matrix.h>
#include <cmn/message.h>
#include <cmn/string.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* allocate new apc molecule struct */
struct apc_mol *mol_apc_new(void) {
  struct apc_mol *p=NULL;
  /* memory allocation */
  p = (struct apc_mol*)malloc(sizeof(struct apc_mol));
  if (!p)
    msg_error("cannot allocate memory for amber prep data struct",1);
  /* initialization */
  p->title = NULL;
  p->file = NULL;
  p->resname = NULL;
  p->n_atoms = 0;
  p->n_loops = 0;
  p->n_imprs = 0;
  p->loop = NULL;
  p->impr = NULL;
  p->atom = NULL;
  return(p);
  }

/* allocate array of apc atom data structures

   n - number of atoms */
struct apc_atom *mol_apc_atom_new(unsigned n) {
  unsigned i,j;
  struct apc_atom *p=NULL;
  /* memory allocation */
  p = (struct apc_atom*)malloc(n*sizeof(struct apc_atom));
  if (!p)
    msg_error("cannot allocate memory for amber prep atom data struct",1);
  /* initialization */
  for (i=0; i<n; i++) {
    p[i].name = NULL;
    p[i].type = NULL;
    p[i].tree = 0;
    p[i].id = 0;
    p[i].charge = 0.0;
    for (j=0; j<3; j++)
      p[i].coord[j] = 0.0;
    }
  return(p);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for apc molecule struct
 
   p - pointer to amber prep struct */
void mol_apc_free(struct apc_mol *p) {
  if (p) {
    p->loop = mat_ufree(p->loop,p->n_loops);
    p->impr = mat_ufree(p->impr,p->n_imprs);
    mol_apc_atom_free(p->atom,p->n_atoms);
    free(p);
    }
  }

/* free memory allocater for array of apc atom data structures

   p - array of atoms
   n - number of atoms */
void mol_apc_atom_free(struct apc_atom *p, unsigned n) {
  unsigned i;
  if (p) {
    for (i=0; i<n; i++) {
      str_free(p[i].name);
      str_free(p[i].type);
      }
    free(p);
    }
  }

/* -------------------------------------------------------------------------- */

/* copy apc atom data

   a0 - the source atom data structure
   a1 - the new atom data structure */
void mol_apc_atom_copy(struct apc_atom *a0, struct apc_atom *a1) {
  unsigned i;
  str_free(a1->name);
  a1->name = str_copy_new(a0->name);
  str_free(a1->type);
  a1->type = str_copy_new(a0->type);
  a1->tree = a0->tree;
  a1->id = a0->id;
  a1->charge = a0->charge;
  for (i=0; i<3; i++)
    a1->coord[i] = a0->coord[i];
  }

/* -------------------------------------------------------------------------- */
