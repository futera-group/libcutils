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

/* return number of atoms in pdb structure

   p - pointer to pdb molecular data struct */
unsigned mol_pdb_atom_n(struct pdb_mol *p) {
  unsigned i,n = 0;
  for (i=0; i<p->n_res; i++)
    n += p->res[i].n_atoms;
  return(n);
  }

/* find residuum and atom ID for given atom 
 
   p - pointer to pdb molecular data struct
   x - sequentional ID of the atom
   r - calculated ID of residuum
   a - calculated ID of atom in the residuum */
short mol_pdb_atom_id(struct pdb_mol *p, unsigned x, unsigned *r, unsigned *a) {
  unsigned i,j,n=0;
  for (i=0; i<p->n_res; i++)
    for (j=0; j<p->res[i].n_atoms; j++)
      if (n++==x) {
        (*r) = i;
        (*a) = j;
        return(1);
        }
  return(0);
  }

/* -------------------------------------------------------------------------- */

/* set unique pdb atomic name using the sequential number
 
   s - storage for the name
   a - atomic number
   n - number of a-atoms in the molecule */
void mol_pdb_atom_name_id(char *s, unsigned a, unsigned n) {
  unsigned m = n;
  if (!n)
    sprintf(s," %s",atom_name(a));
  else if (n<10)
    sprintf(s," %s%d",atom_name(a),n);
  else {
    while (m>=10)
      m/=10;
    sprintf(s,"%d%s%d",m,atom_name(a),(n>=100 ? n%100 : n%10));
    }
  s[4]='\0';
  }

/* set pdb atomic name and fix the letter position

   p    - the PDB molecular data
   ir   - residuum ID
   ia   - atom ID
   name - name of the atom */
void mol_pdb_atom_name_pos(struct pdb_mol *p, unsigned ir,
  unsigned ia, char *name) {
  char *s;
  s = atom_name_to_pdb(name);
  str_free(p->res[ir].atom[ia].name);
  p->res[ir].atom[ia].name = str_copy_new(s);
  }

/* -------------------------------------------------------------------------- */
