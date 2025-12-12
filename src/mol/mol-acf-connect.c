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

#include <cmn/matrix.h>
#include <cmn/vector.h>
#include "mol/distance.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* find minimal interatomic distance between specified fragments
 
   c     - pointer to acf molecular data struct
   a     - array with fragmentation IDs
   f1,f2 - fragment IDs
   a1,a2 - atom IDs (output) */
double mol_acf_connect_fdmin(struct acf_mol *c, unsigned *a, unsigned f1,
  unsigned f2, unsigned *a1, unsigned *a2) {
  double d,dm = 9.0E+90;
  unsigned i,j;
  for (i=0; i<c->n_atoms; i++)
    if (a[i]==f1)
      for (j=0; j<c->n_atoms; j++)
        if (a[j]==f2) {
          d = dist_r(c->atom[i].coord,c->atom[j].coord,3);
          if (d<dm) {
            (*a1) = i;
            (*a2) = j;
            dm = d;
            }
          }
  return(dm);
  }

/* define additional bonds to connect fragments of atoms
 
   c - pointer to acf molecular data struct
   a - array with fragmentation IDs
   n - number of fragments  */
void mol_acf_connect_bonds(struct acf_mol *c, unsigned *a, unsigned n) {
  unsigned i,j,k,mn_f1,mn_f2,mn_a1,mn_a2,a1,a2;
  double mn_d,d;
  while (n>1) {
    mn_d = 9.0E+90;
    mn_f1 = mn_f2 = mn_a1 = mn_a2 = 0;
    /* find minimal distance */
    for (i=0; i<n; i++)
      for (j=i+1; j<n; j++) {
        d = mol_acf_connect_fdmin(c,a,i,j,&a1,&a2);
        if (d<mn_d) {
          mn_f1 = i;
          mn_a1 = a1;
          mn_f2 = j;
          mn_a2 = a2;
          mn_d = d;
          }
        }
    /* define new bond */
    c->bond = mat_uresize_r(c->bond,2,c->n_bonds,c->n_bonds+1);
    c->bond[c->n_bonds][0] = mn_a1;
    c->bond[c->n_bonds][1] = mn_a2;
    c->n_bonds++;
    /* merge fragment IDs */
    for (k=0; k<c->n_atoms; k++)
      if (a[k]==mn_f2)
        a[k] = mn_f1;
    if (mn_f2<(n-1))
      for (k=0; k<c->n_atoms; k++)
        if (a[k]==(n-1))
          a[k] = mn_f2;
    n--;
    }
  }

/* check fragmentation and define additional bonds to connect all atoms
 
   c - pointer to acf molecular data struct */
void mol_acf_connect(struct acf_mol *c) {
  unsigned i,j,id = 1,rr,*at,*a2;
  /* auxiliary array */
  at = vec_ualloc(c->n_atoms);
  vec_uset(at,0,c->n_atoms);
  /* color fragments */
  for (i=0; i<c->n_bonds; i++) {
    if (at[c->bond[i][0]]==0 && at[c->bond[i][1]]==0) {
      at[c->bond[i][0]] = id;
      at[c->bond[i][1]] = id;
      id++;
      }
    else if (at[c->bond[i][0]]==0)
      at[c->bond[i][0]] = at[c->bond[i][1]];
    else if (at[c->bond[i][1]]==0)
      at[c->bond[i][1]] = at[c->bond[i][0]];
    else {
      rr = at[c->bond[i][0]];
      for (j=0; j<c->n_atoms; j++)
        if (at[j]==rr)
          at[j] = at[c->bond[i][1]];
      }
    }
  /* count fragments */
  id = 1;
  a2 = vec_ucopy_new(at,c->n_atoms);
  vec_usort(a2,c->n_atoms); 
  for (i=1; i<c->n_atoms; i++)
    if (a2[i-1]!=a2[i]) {
      for (j=0; j<c->n_atoms; j++)
        if (at[j]==a2[i-1])
          at[j] = id-1;
      id++;
      }
  for (i=0; i<c->n_atoms; i++)
    if (at[i]==a2[c->n_atoms-1])
      at[i] = id-1;
  /* define new bonds */
  if (id>1)
   mol_acf_connect_bonds(c,at,id);
  /* clean memory */
  vec_ufree(at);
  vec_ufree(a2);
  }

/* -------------------------------------------------------------------------- */
