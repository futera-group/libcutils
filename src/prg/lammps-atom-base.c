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

/* allocate memory for LAMMPS atomic data
   
   n - number of the atoms */
struct lammps_atom *lammps_atom_new(unsigned n) {
  unsigned i;
  struct lammps_atom *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct lammps_atom));
  if (!d) 
    msg_error("cannot allocate memory for LAMMPS atomic data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].type = 0;
    d[i].mol = 0;
    d[i].charge = 0.0;
    d[i].crd = NULL;
    d[i].vel = NULL;
    d[i].rvl = NULL;
    d[i].avl = NULL;
    d[i].ang = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for LAMMPS atomic data 
 
   d - the atomic data
   n - number of the atoms */
struct lammps_atom *lammps_atom_free(struct lammps_atom *d, unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++) {
      d[i].crd = vec_ffree(d[i].crd);
      d[i].vel = vec_ffree(d[i].vel);
      d[i].rvl = vec_ffree(d[i].rvl);
      d[i].avl = vec_ffree(d[i].avl);
      d[i].ang = vec_ffree(d[i].ang);
      }
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
