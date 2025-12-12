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
#include "mol/cell.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* allocate new general molecule struct */
struct gen_mol *mol_gen_new(void) {
  struct gen_mol *m = NULL;
  /* allocate memory */
  m = (struct gen_mol*)malloc(sizeof(struct gen_mol));
  if (!m)
    msg_error("cannot allocate memory for general molecule data struct",1);
  /* initialization */
  m->title = NULL;
  m->n_atoms = 0;
  m->atom = NULL;
  m->cell = NULL;
  return(m);
  }

/* allocate array of general atom data structs

   n - number of atoms */
struct gen_atom* mol_gen_atom_new(unsigned n) {
  unsigned i,j;
  struct gen_atom *a = NULL;
  /* allocate memory */
  a = (struct gen_atom*)malloc(n*sizeof(struct gen_atom));
  if (!a)
    msg_error("cannot allocate memory of general atom data array",1);
  /* initialization */
  for (i=0; i<n; i++) {
    a[i].num = 0;
    a[i].charge = 0.0;
    for (j=0; j<3; j++)
      a[i].coord[j] = 0.0;
    }
  return(a);
  }

/* -------------------------------------------------------------------------- */

/* free memory allocated for general molecule struct

   m - pointer to general molecular data struct */
void mol_gen_free(struct gen_mol *m) {
  if (m) {
    str_free(m->title);
    mol_gen_atom_free(m->atom);
    cell_free(m->cell);
    free(m);
    }
  }

/* free memory allocated for general atom data 

   a - atom data */
void mol_gen_atom_free(struct gen_atom *a) {
  if (a)
    free(a);
  }

/* -------------------------------------------------------------------------- */
