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
#include "mol/atom.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* convert apc molecular data struct to acf struct
 
   p - pointer to amber prep data struct */
struct acf_mol *mol_apc_acf(struct apc_mol *p) {
  char sym[5];
  unsigned i,j;
  struct acf_mol *m;
  /* apc structure */
  m = mol_acf_new();
  m->n_atoms = p->n_atoms;
  m->atom = mol_acf_atom_new(m->n_atoms);
  /* convert molecular representation */
  m->name = str_copy_new(p->resname);
  for (i=0; i<m->n_atoms; i++) {
    sprintf(sym," %-3s",p->atom[i].name);
    m->atom[i].name = str_copy_new(sym);
    m->atom[i].type = str_copy_new(p->atom[i].type);
    for (j=0; j<3; j++)
      m->atom[i].coord[j] = p->atom[i].coord[j];
    m->atom[i].num = atom_num_pdb(p->atom[i].name);
    m->atom[i].charge = p->atom[i].charge;
    m->charge += m->atom[i].charge;
    }
  /* chemical formula */
  mol_acf_set_formula(m);
  /* interatomic bonds */
  mol_acf_set_bonds(m);
  /* atomic types */
  mol_acf_type_gaff(m);
  return(m);
  }

/* -------------------------------------------------------------------------- */
