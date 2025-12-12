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
#include "prg/molden.h"

/* -------------------------------------------------------------------------- */

/* allocate new molden data struct */
struct mld_dat *mld_dat_new(void) {
  struct mld_dat *g = NULL;
  /* memory allocation */
  g = (struct mld_dat*)malloc(sizeof(struct mld_dat));
  if (!g)
    msg_error("cannot allocate new molden data struct",1);
  /* initialization */
  g->n_atoms = 0;
  g->n_vib_modes = 0;
  g->atom = NULL;
  g->freq = NULL;
  return(g);
  }

/* free memory allocated for molden data struct
   
   g - pointer to molden data struct */
struct mld_dat* mld_dat_free(struct mld_dat *g) {
  if (g) {
    g->atom = mld_atom_free(g->atom,g->n_atoms);
    g->freq = mld_freq_free(g->freq,g->n_vib_modes,g->n_atoms);
    free(g);
    }
  return(NULL);
  }

/* -------------------------------------------------------------------------- */
