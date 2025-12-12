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
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* allocate vector of atom gaussian data struct
 
   n - number of atoms */
struct gauss_atom *gauss_atom_new(unsigned n) {
  unsigned i,j;
  struct gauss_atom *a = NULL;
  /* memory allocation */
  a = (struct gauss_atom*)malloc(n*sizeof(struct gauss_atom));
  if (!a)
    msg_error("cannot allocate memory for gaussian atom data structure",1);
  /* initilalization */
  for (i=0; i<n; i++) {
    for (j=0; j<3; j++) {
      a[i].coord[j] = 0.0;
      a[i].coord_fix[j] = 0.0;
      a[i].grad[j] = 0.0;
      }
    for (j=0; j<5; j++)
      a[i].charge[j] = 0.0;
    a[i].num = 0;
    a[i].mass = 0.0;
    a[i].type_name[0][0] = '\0';
    a[i].type_name[1][0] = '\0';
    a[i].type_id[0] = 0;
    a[i].type_id[1] = 0;
    a[i].weight = 0;
    a[i].res_info = 0;
    a[i].res_num = 0;
    a[i].frag_info = 0;
    a[i].nuc_spin = 0;
    a[i].nuc_qmom = 0.0;
    a[i].nuc_gfac = 0.0;
    }
  return(a);
  }

/* free memory allocated for gaussian atom struct 

   g - pointer to gaussian atom data struct
   n - number of atoms */
struct gauss_atom* gauss_atom_free(struct gauss_atom *g, unsigned n) {
  if (g)
    free(g);
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
