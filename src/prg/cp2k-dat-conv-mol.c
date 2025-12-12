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

#include <cmn/string.h>
#include <mol/cell.h>
#include <mol/molec.h>
#include "prg/cp2k.h"

/* -------------------------------------------------------------------------- */

/* convert cp2k structure to generic molecular format
 
   d - cp2k data struct */
struct gen_mol *cp2k_dat_conv_mol(struct cp2k_dat *d) {
  unsigned i,j;
  struct gen_mol *m = NULL;
  m = mol_gen_new();
  /* cell */
  m->cell = cell_copy_new(d->cell);
  /* structure */
  m->n_atoms = d->n_atoms;
  m->atom = mol_gen_atom_new(m->n_atoms);
  m->title = str_copy_new("CP2K Input Structure");
  for (i=0; i<d->n_atoms; i++) {
    for (j=0; j<3; j++)
      m->atom[i].coord[j] = d->atom[i].crd[j];
    m->atom[i].num = d->atom[i].num;
    }
  return(m);
  }

/* -------------------------------------------------------------------------- */
