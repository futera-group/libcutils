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

/* allocate memory for DL_POLY RDF parameters
   
   n - number of the RDF pairs */
struct dlpoly_prm_rdf *dlpoly_prm_rdf_new(unsigned n) {
  unsigned i;
  struct dlpoly_prm_rdf *d = NULL;
  /* memory allocation */
  d = malloc(n*sizeof(struct dlpoly_prm_rdf));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY RDF pair data",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].at1 = NULL;
    d[i].at2 = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for DL_POLY RDF pair data 
 
   d - the RDF pair data
   n - number of the data structures */
struct dlpoly_prm_rdf *dlpoly_prm_rdf_free(struct dlpoly_prm_rdf *d,
  unsigned n) {
  unsigned i;
  if (d) {
    for (i=0; i<n; i++) {
      d[i].at1 = str_free(d[i].at1);
      d[i].at2 = str_free(d[i].at2);
      }
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
