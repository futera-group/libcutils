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
#include <cmn/matrix.h>
#include <cmn/message.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* allocate vector of structure gaussian data struct
 
   n - number of structures */
struct gauss_geom *gauss_geom_new(unsigned n) {
  unsigned i;
  struct gauss_geom *g = NULL;
  /* memory allocation */
  g = (struct gauss_geom*)malloc(n*sizeof(struct gauss_geom));
  if (!g)
    msg_error("cannot allocate memory for gaussian geom data struct",1);
  /* initialization */
  for (i=0; i<n; i++) {
    g[i].energy = 0.0;
    g[i].coord = NULL;
    }
  return(g);
  }

/* free memory allocated for gaussian geom data struct

   g  - pointer to vector of gaussian geometry data structs 
   ng - number of structures
   na - number of atoms in each structure */
struct gauss_geom* gauss_geom_free(struct gauss_geom *g, unsigned ng,
  unsigned na) {
  unsigned i;
  if (g) {
    for (i=0; i<ng; i++)
      g[i].coord = mat_ffree(g[i].coord,na);
    free(g);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
