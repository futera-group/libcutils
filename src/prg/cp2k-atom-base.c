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
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* allocate and initialize array of atoms
 
   n - number of the atoms */
struct cp2k_atom *cp2k_atom_new(unsigned n) {
  unsigned i,j;
  struct cp2k_atom *d = NULL;
  /* memory allocation */
  d = (struct cp2k_atom*)malloc(n*sizeof(struct cp2k_atom));
  if (!d)
    msg_error("cannot allocate memory for cp2k atom array",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].num = 0;
    d[i].kind = 0;
    for (j=0; j<3; j++)
      d[i].charge[j] = 0.0;
    d[i].mass = 0.0;
    for (j=0; j<3; j++)
      d[i].crd[j] = 0.0;
    }
  return(d);
  }

/* free memory allocated for cp2k atom array
 
   d - the atom array */
struct cp2k_atom* cp2k_atom_free(struct cp2k_atom *d) {
  if (d) {
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
