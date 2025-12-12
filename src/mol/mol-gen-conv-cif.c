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
#include <cmn/vector.h>
#include "mol/atom.h"
#include "mol/cell.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert generic molecular data struct to cif format

   m - pointer to generic molecular data struct */
struct cif_mol* mol_gen_cif(struct gen_mol *m) {
  char s[1024];
  unsigned i,j,a[120];
  struct cif_mol *c;
  /* initialization */
  vec_uset(a,1,120);
  /* cif structure */
  c = mol_cif_new();
  c->cell = cell_copy_new(m->cell);
  c->n_atoms = m->n_atoms;
  c->atom = mol_cif_atom_new(c->n_atoms);
  for (i=0; i<c->n_atoms; i++) {
    c->atom[i].num = m->atom[i].num;
    sprintf(s,"%s%d",atom_name(c->atom[i].num),a[c->atom[i].num]++);
    c->atom[i].name = str_copy_new(s);
    for (j=0; j<3; j++)
      c->atom[i].coord[j] = m->atom[i].coord[j];
    }
  return(c);
  }

/* -------------------------------------------------------------------------- */
