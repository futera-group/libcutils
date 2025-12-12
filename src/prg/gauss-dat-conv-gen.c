/******************************************************************************\
 *                                                                            * 
 *  Libcutils - library of C function                                         * 
 *                                                                            *
 *  Version:             3.4                                                  * 
 *  Date:                06/03/2020                                           *
 *                                                                            * 
 *  Author:              Zdenek Futera                                        * 
 *                                                                            * 
 *  Address:             University of South Bohemia                          * 
 *                       Faculty of Science, Institute of Physics             * 
 *                       Branisovska 1760, 370 05 Ceske Budejovice            * 
 *                       Czech Republic                                       * 
 *                                                                            * 
 *  E-Mail:              zfutera@prf.jcu.cz                                   * 
 *                                                                            * 
\******************************************************************************/

#include <cmn/string.h>
#include <cmn/units.h>
#include <mol/molec.h>
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* convert gaussian data to generic molecular file format
 
   g - the gaussian data */
struct gen_mol* gauss_dat_gen(struct gauss_dat *g) {
  unsigned i,j;
  struct gen_mol *d;
  d = mol_gen_new();
  d->n_atoms = g->n_atoms;
  d->atom = mol_gen_atom_new(d->n_atoms);
  d->title = str_copy_new(g->job_title);
  for (i=0; i<d->n_atoms; i++) {
    d->atom[i].num = g->atom[i].num;
    for (j=0; j<3; j++)
      d->atom[i].coord[j] = g->atom[i].coord[j]*CONV_B_ANG;
    d->atom[i].charge = g->atom[i].charge[GAUSS_CH_MULL];
    }
  return(d);
  }

/* -------------------------------------------------------------------------- */
