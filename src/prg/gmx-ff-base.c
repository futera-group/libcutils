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
#include <cmn/vector.h>
#include "prg/gromacs.h"

/* -------------------------------------------------------------------------- */

/* Allocate memory for gromacs force-field data structure */
struct gmx_ff* gmx_ff_new(void) {
  struct gmx_ff *d = NULL;
  /* memory allocation */
  d = (struct gmx_ff*)malloc(sizeof(struct gmx_ff));
  if (!d)
    msg_error("cannot allocate memory for force-field data",1);
  /* initialization */
  d->file = NULL;
  d->n_defs = 0;
  d->def = NULL;
  d->nb_func = 0;
  d->comb_rule = 0;
  d->gen_pairs = 0;
  d->scale_lj = 0.0;
  d->scale_qq = 0.0;
  d->dat = gmx_top_new();
  return(d);
  }

/* Clean memory allocated in gromacs force-field data structure

   d - the data structure */
void gmx_ff_clean(struct gmx_ff *d) {
  d->file = str_free(d->file);
  vec_sfree(d->def,d->n_defs);
  d->n_defs = 0;
  d->nb_func = 0;
  d->comb_rule = 0;
  d->gen_pairs = 0;
  d->scale_lj = 0.0;
  d->scale_qq = 0.0;
  d->dat = gmx_top_free(d->dat);
  }

/* Clean memory allocated for gromacs force-field data structure

   d - the data array */
struct gmx_ff* gmx_ff_free(struct gmx_ff *d) {
  if (d) {
    gmx_ff_clean(d);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* Copy gromacs data from one force-field data structure to another

   d1 - destination data structure
   d2 - source data structure */
void gmx_ff_copy(struct gmx_ff *d1, struct gmx_ff *d2) {
  unsigned i;
  if (d1 && d2) {
    d1->file = str_copy_new(d2->file);
    d1->n_defs = d2->n_defs;
    d1->def = vec_talloc(sizeof(char*),d1->n_defs);
    for (i=0; i<d1->n_defs; i++)
      d1->def[i] = str_copy_new(d2->def[i]);
    d1->nb_func = d2->nb_func;
    d1->comb_rule = d2->comb_rule;
    d1->gen_pairs = d2->gen_pairs;
    d1->scale_lj = d2->scale_lj;
    d1->scale_qq = d2->scale_qq;
    d1->dat = gmx_top_copy_new(d2->dat);
    }
  }

/* Create copy of gromacs force-field data structure

   d - the data array
   n - length of the array */
struct gmx_ff* gmx_ff_copy_new(struct gmx_ff *d) {
  struct gmx_ff *t = NULL;
  if (d) {
    t = gmx_ff_new();
    gmx_ff_copy(t,d);
    }
  return(t);
  }

/* -------------------------------------------------------------------------- */
