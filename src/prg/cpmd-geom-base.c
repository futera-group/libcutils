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
#include "prg/cpmd.h"

/* -------------------------------------------------------------------------- */

/* allocate array of geom data structs

   n - number of structures
   m - number of coordinates */
struct cpmd_geom *cpmd_geom_new(unsigned n, unsigned m) {
  unsigned i;
  struct cpmd_geom *d = NULL;
  /* memory allocation */
  d = (struct cpmd_geom*)malloc(n*sizeof(struct cpmd_geom));
  if (!d)
    msg_error("cannot allocate memory for cpmd geometry data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].energy = 0.0;
    d[i].coord = mat_falloc(m,3);
    mat_fset(d[i].coord,0.0,m,3);
    }
  return(d);
  }

/* free memory allocated for array of geom data structs
 
   d - pointer to the array
   n - number of the structs
   m - number of coordinates */
struct cpmd_geom *cpmd_geom_free(struct cpmd_geom *d, unsigned n, unsigned m) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++)
      d[i].coord = mat_ffree(d[i].coord,m);
    free(d);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
