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
#include <cmn/string.h>
#include "prg/dlpoly.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for DL_POLY atomic data
   
   n - number of the atoms */
struct dlpoly_atom *dlpoly_atom_new(unsigned n) {
  unsigned i;
  struct dlpoly_atom *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct dlpoly_atom));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY atomic data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].name = NULL;
    d[i].n_atoms = 0;
    d[i].frozen = 0;
    d[i].mass = 0.0;
    d[i].charge = 0.0;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for DL_POLY atomic data 
 
   d - the atomic data
   n - number of atoms */
struct dlpoly_atom *dlpoly_atom_free(struct dlpoly_atom *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++)
      d[i].name = str_free(d[i].name);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
