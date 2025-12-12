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
#include <cmn/vector.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* allocate memory for array electronic-coupling data blocks

   n - number of the blocks */
struct cp2k_ecpl_block *cp2k_ecpl_block_new(unsigned n) {
  unsigned i;
  struct cp2k_ecpl_block *d = NULL;
  /* memory allocation */
  d = (struct cp2k_ecpl_block*)malloc(n*sizeof(struct cp2k_ecpl_block));
  if (!d)
    msg_error("cannot allocate memory for cp2k coupling blocks",1);
  /* initialization */
  for (i=0; i<n; i++) {
    d[i].n_atoms = 0;
    d[i].atom_id = NULL;
    d[i].n_states = NULL;
    d[i].energy = NULL;
    d[i].ao_pj = NULL;
    d[i].coupling = NULL;
    }
  return(d);
  }

/* free memory allocated for array of cp2k electronic-coupling data blocks
 
   d        - the electronic-coupling blocks
   n_blocks - number of the blocks
   n_spins  - number of spin components */
struct cp2k_ecpl_block* cp2k_ecpl_block_free(struct cp2k_ecpl_block *d,
  unsigned n_blocks, unsigned n_spins) {
  unsigned ib1,ib2,is1,ix,id;
  if (d) {
    for (ib1=0; ib1<n_blocks; ib1++) {
      d[ib1].atom_id = vec_ufree(d[ib1].atom_id);
      d[ib1].energy = mat_ffree(d[ib1].energy,n_spins);
      if (d[ib1].ao_pj) {
        for (ix=0; ix<n_spins; ix++)
          d[ib1].ao_pj[ix] = mat_ffree(d[ib1].ao_pj[ix],d[ib1].n_states[ix]);
        d[ib1].ao_pj = vec_tfree(d[ib1].ao_pj);
        }
      if (d[ib1].coupling) {
        for (ix=0; ix<n_spins; ix++) {
          for (ib2=ib1+1; ib2<n_blocks; ib2++) {
            id = ib2-ib1-1;
            for (is1=0; is1<d[ib1].n_states[ix]; is1++)
              d[ib1].coupling[ix][id][is1] = 
                vec_ffree(d[ib1].coupling[ix][id][is1]);
            d[ib1].coupling[ix][id] = vec_tfree(d[ib1].coupling[ix][id]);
            }
          d[ib1].coupling[ix] = vec_tfree(d[ib1].coupling[ix]);
          }
        }
      d[ib1].n_states = vec_ufree(d[ib1].n_states);
      d[ib1].coupling = vec_tfree(d[ib1].coupling);
      }
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */

/* allocate memory for electronic-coupling data structure */
struct cp2k_ecpl *cp2k_ecpl_new(void) {
  struct cp2k_ecpl *d;
  /* memory allocation */
  d = (struct cp2k_ecpl*)malloc(sizeof(struct cp2k_ecpl));
  if (!d)
    msg_error("cannot allocate memory for cp2k couplings",1);
  /* initialization */
  d->n_blocks = 0;
  d->n_states = NULL;
  d->energy = NULL;
  d->block = NULL;
  return(d);
  }

/* free memory allocated for cp2k electronic-coupling data structure
 
   d       - the electronic coupling data
   n_spins - number of spin components */
struct cp2k_ecpl* cp2k_ecpl_free(struct cp2k_ecpl *d, unsigned n_spins) {
  if (d) {
    d->n_states = vec_ufree(d->n_states);
    d->energy = mat_ffree(d->energy,n_spins);
    d->block = cp2k_ecpl_block_free(d->block,d->n_blocks,n_spins);
    free(d);
    d = NULL;
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
