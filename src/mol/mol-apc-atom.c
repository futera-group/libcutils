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

#include <cmn/message.h>
#include <cmn/string.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* return ID of the atom specified by name
 
   p   - the amber prep data structure
   sym - name of the atom */
unsigned mol_apc_atom_id(struct apc_mol *p, char *sym) {
  unsigned i;
  for (i=0; i<p->n_atoms; i++)
    if (str_compare(sym,p->atom[i].name))
      return(i);
  msg_error_f("undefined atom name \"%s\" specified for ID search",1,sym);
  return(0);
  }

/* return ID of specified atom
 
   p - the amber prep data structure */
unsigned mol_apc_atom_id2id(struct apc_mol *p, unsigned a) {
  unsigned i;
  for (i=0; i<p->n_atoms; i++)
    if (p->atom[i].id==a)
      return(i);
  return(0);
  }

/* reset all atom IDs as well as loop and improper pointers
 
   m - the amber prep data structure */
void mol_apc_atom_id_set(struct apc_mol *m) {
  unsigned i,j;
  /* renumber loops */
  for (i=0; i<m->n_loops; i++)
    for (j=0; j<2; j++)
      m->loop[i][j] = mol_apc_atom_id2id(m,m->loop[i][j]);
  /* renumber impropers */
  for (i=0; i<m->n_imprs; i++)
    for (j=0; j<4; j++)
      m->impr[i][j] = mol_apc_atom_id2id(m,m->loop[i][j]);
  /* reset IDs */
  for (i=0; i<m->n_atoms; i++)
    m->atom[i].id = i+4;
  }

/* -------------------------------------------------------------------------- */
