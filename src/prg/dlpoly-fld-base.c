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

/* allocate memory for DL_POLY force-field parameters */
struct dlpoly_fld *dlpoly_fld_new(void) {
  struct dlpoly_fld *d = NULL;
  /* memory allocation */
  d = malloc(sizeof(struct dlpoly_fld));
  if (!d) 
    msg_error("cannot allocate memory for DL_POLY force-field data",1);
  /* initialization */
  d->header = NULL;
  d->units = 0;
  d->n_molecules = 0;
  d->mol = NULL;
  d->n_vdws = 0;
  d->vdw = NULL;
  d->n_metals = 0;
  d->metal = NULL;
  d->n_rdfs = 0;
  d->rdf = NULL;
  d->n_tersoffs = 0;
  d->ters = NULL;
  d->n_tbps = 0;
  d->tbp = NULL;
  d->n_fbps = 0;
  d->fbp = NULL;
  d->n_ext_fields = 0;
  d->efld = NULL;
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for DL_POLY force-field potential data 
 
   d - the force-field data */
struct dlpoly_fld *dlpoly_fld_free(struct dlpoly_fld *d) {
  if (d) {
    d->header = str_free(d->header);
    d->mol = dlpoly_mol_free(d->mol,d->n_molecules);
    d->vdw = dlpoly_prm_vdw_free(d->vdw,d->n_vdws);
    d->metal = dlpoly_prm_metal_free(d->metal,d->n_metals);
    d->rdf = dlpoly_prm_rdf_free(d->rdf,d->n_rdfs);
    d->ters = dlpoly_prm_ters_free(d->ters,d->n_tersoffs);
    d->tbp = dlpoly_prm_tbp_free(d->tbp,d->n_tbps);
    d->fbp = dlpoly_prm_fbp_free(d->fbp,d->n_fbps);
    d->efld = dlpoly_prm_efld_free(d->efld);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
