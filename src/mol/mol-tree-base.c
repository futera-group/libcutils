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
#include <cmn/message.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for new tree atomic data structure */
struct mol_tree_atom* mol_tree_atom_new(void) {
  struct mol_tree_atom *t = NULL;
  /* memory allocation */
  t = (struct mol_tree_atom*)malloc(sizeof(struct mol_tree_atom));
  if (!t) 
    msg_error("cannot allocate memory for new tree atomic data structure",1);
  /* initialization */
  t->id = 0;
  t->tree = 0;
  return(t);
  }

/* free memory allocated for a tree atomic data structure

   t - the tree atomic data structure */
void mol_tree_atom_free(struct mol_tree_atom *t) {
  if (t)
    free(t);
  }

/* -------------------------------------------------------------------------- */
