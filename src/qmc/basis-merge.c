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
#include "qmc/basis.h"

/* -------------------------------------------------------------------------- */

/* merge two basis sets to one basis set data struct
 
   b1,b2 - the basis sets for merging */
struct basis* basis_merge(struct basis *b1, struct basis *b2) {
  unsigned i;
  struct basis *b;
  b = basis_new();
  /* merge centers */
  b->n_centers = (b1->n_centers+b2->n_centers);
  b->center = basis_center_new(b->n_centers);
  for (i=0; i<b1->n_centers; i++)
    basis_center_copy(b->center+i,b1->center+i);
  for (i=0; i<b2->n_centers; i++)
    basis_center_copy(b->center+b1->n_centers+i,b2->center+i);
  /* set global parameters */
  basis_def_parm(b);
  b->n_ibfce = (b1->n_ibfce+b2->n_ibfce);
  /* check consistency */
  if (b->n_bfce!=(b1->n_bfce+b2->n_bfce))
    msg_error("logic error in merged basis set (n_bfce)",1);
  return(b);
  }

/* -------------------------------------------------------------------------- */
