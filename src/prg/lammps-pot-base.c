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
#include <cmn/vector.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for LAMMPS potential data
   
   n - number of types */
struct lammps_pot *lammps_pot_new(unsigned n) {
  unsigned i;
  struct lammps_pot *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct lammps_pot));
  if (!d) 
    msg_error("cannot allocate memory for LAMMPS potential data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].type = 0;
    d[i].style = 0;
    d[i].n_ids = 0;
    d[i].id = NULL;
    d[i].n_coeffs = 0;
    d[i].coeff = NULL;
    }
  return(d);
  }

/* free memory allocated for LAMMPS potential data 
 
   d - potential data
   n - number of the terms */
struct lammps_pot *lammps_pot_free(struct lammps_pot *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++) {
      d[i].id = vec_ufree(d[i].id);
      d[i].coeff = vec_ffree(d[i].coeff);
      }
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
