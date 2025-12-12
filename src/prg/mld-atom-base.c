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
#include "prg/molden.h"

/* -------------------------------------------------------------------------- */

/* allocate vector of atom molden data struct
 
   n - number of atoms */
struct mld_atom *mld_atom_new(unsigned n) {
  unsigned i,j;
  struct mld_atom *a = NULL;
  /* memory allocation */
  a = (struct mld_atom*)malloc(n*sizeof(struct mld_atom));
  if (!a)
    msg_error("cannot allocate memory for molden atom data structure",1);
  /* initilalization */
  for (i=0; i<n; i++) {
    for (j=0; j<3; j++)
      a[i].crd[j] = 0.0;
    a[i].sym = NULL;
    a[i].num = 0;
    }
  return(a);
  }

/* free memory allocated for molden atom struct 

   m - pointer to molden atom data struct
   n - number of atoms */
struct mld_atom* mld_atom_free(struct mld_atom *m, unsigned n) {
  if (m)
    free(m);
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
