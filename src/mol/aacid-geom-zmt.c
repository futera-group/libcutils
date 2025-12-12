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
#include <cmn/list.h>
#include "mol/aacid.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* Create structure of a single amino acid (Z-matrix format)
 
   id   - internal ID of the amino acid
   form - protonation form of the amino acid */
struct zmt_mol* aacid_geom_zmt(short id, short form) {
  unsigned i;
  struct list *z;
  struct ldata *p;
  struct zmt_mol *t;
  /* define amino-acid structure */
  z = aacid_geom_def(id,form);
  /* create Z-matrix */
  t = mol_zmt_new();
  t->n_atoms = z->num;
  t->atom = mol_zmt_atom_new(t->n_atoms);
  for (p=z->first,i=0; p; p=p->l_next,i++)
    mol_zmt_atom_copy((struct zmt_atom*)p->l_data,t->atom+i);
  /* clean memory */
  list_free(z,free);
  return(t);
  }

/* -------------------------------------------------------------------------- */
