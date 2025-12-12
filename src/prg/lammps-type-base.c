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
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for LAMMPS particle type data
   
   n - number of atomic types */
struct lammps_type *lammps_type_new(unsigned n) {
  unsigned i;
  struct lammps_type *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct lammps_type));
  if (!d) 
    msg_error("cannot allocate memory for LAMMPS particle type data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].name = NULL;
    d[i].num = 0;
    d[i].mass = 0.0;
    }
  return(d);
  }

/* free memory allocated for LAMMPS particle type data 
 
   d - particle type data
   n - number of types */
struct lammps_type *lammps_type_free(struct lammps_type *d, unsigned n) {
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
