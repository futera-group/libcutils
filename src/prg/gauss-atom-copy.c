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

#include <string.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* copy data from one gaussian data struct to another 
 
   a1 - the destination data struct where the data are writen
   a2 - the source data struct where the data are read */
void gauss_atom_copy(struct gauss_atom *a1, struct gauss_atom *a2) {
  unsigned i;
  for (i=0; i<3; i++) {
    a1->coord[i] = a2->coord[i];
    a1->coord_fix[i] = a2->coord_fix[i];
    a1->grad[i] = a2->grad[i];
    }
  a1->num = a2->num;
  for (i=0; i<5; i++)
    a1->charge[i] = a2->charge[i];
  a1->mass = a2->mass;
  for (i=0; i<2; i++) {
    strncpy(a1->type_name[i],a2->type_name[i],25);
    a1->type_id[i] = a2->type_id[i];
    }
  a1->weight = a2->weight;
  a1->res_info = a2->res_info;
  a1->res_num = a2->res_num;
  a1->frag_info = a2->frag_info;
  a1->nuc_spin = a2->nuc_spin;
  a1->nuc_qmom = a2->nuc_qmom;
  a1->nuc_gfac = a2->nuc_gfac;
  }

/* -------------------------------------------------------------------------- */
