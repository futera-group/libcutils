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
#include <cmn/vector.h>
#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

/* create bond array with atom indices for given atom
  
   t - pointer to amber topology struct
   a - atom ID
   n - number of bonds (output)
   q - vector for equilibrium values (ignored if NULL)
   h - type of bonds */
unsigned *amber_top_bonds_atom(struct amber_top *t, unsigned a, unsigned *n, 
  double **q, short h) {
  unsigned i,*b,nh,na,id1,id2,n_bonds = 0;
  double *d;
  nh = t->pointers[AMBER_POINTER_NBONH];
  na = t->pointers[AMBER_POINTER_NBONA];
  /* bond array */
  b = vec_ualloc(nh+na);
  d = vec_falloc(nh+na);
  /* bonds including H */
  if (h==AMBER_BOND_ALL || h==AMBER_BOND_WITH_H)
    for (i=0; i<nh; i++) {
      id1 = abs(t->bond_h[3*i]/3);
      id2 = abs(t->bond_h[3*i+1]/3);
      if (id1==a) {
        b[n_bonds] = id2;
        d[n_bonds] = t->bond_eq_val[t->bond_h[3*i+2]-1];
        n_bonds++;
        }
      else if (id2==a) {
        b[n_bonds] = id1;
        d[n_bonds] = t->bond_eq_val[t->bond_h[3*i+2]-1];
        n_bonds++;
        }
      }
  /* bonds without H */
  if (h==AMBER_BOND_ALL || h==AMBER_BOND_WITHOUT_H)
    for (i=0; i<na; i++) {
      id1 = abs(t->bond_ah[3*i]/3);
      id2 = abs(t->bond_ah[3*i+1]/3);
      if (id1==a) {
        b[n_bonds] = id2;
        d[n_bonds] = t->bond_eq_val[t->bond_ah[3*i+2]-1];
        n_bonds++;
        }
      else if (id2==a) {
        b[n_bonds] = id1;
        d[n_bonds] = t->bond_eq_val[t->bond_ah[3*i+2]-1];
        n_bonds++;
        }
      }
  /* final size correction */
  if (n_bonds<(na+nh)) {
    b = vec_uresize(b,na+nh,n_bonds);
    d = vec_fresize(d,na+nh,n_bonds);
    }
  if (q)
    (*q)=d;
  else
    vec_ffree(d);
  (*n) = n_bonds;
  return(b);
  }

/* create bond array with atom indices for given residuum
  
   t - pointer to amber topology struct
   r - residuum ID
   n - number of bonds (output)
   q - vector for equilibrium values (ignored if NULL)
   h - type of bonds */
unsigned **amber_top_bonds_res(struct amber_top *t, unsigned r, unsigned *n,
  double **q, short h) {
  unsigned i,**b,nh,na,a0,a1,id1,id2,n_atoms,n_bonds = 0;
  double *d;
  /* number of atoms and bonds */
  n_atoms = (r==(t->pointers[AMBER_POINTER_NRES]-1) ?
    t->pointers[AMBER_POINTER_NATOM]-t->res_pointer[r]+1 :
    t->res_pointer[r+1]-t->res_pointer[r]);
  nh = t->pointers[AMBER_POINTER_NBONH];
  na = t->pointers[AMBER_POINTER_NBONA];
  /* first and last atom index */
  a0 = t->res_pointer[r]-1;
  a1 = a0+n_atoms;
  /* bond array */
  b = mat_ualloc(nh+na,2);
  d = vec_falloc(nh+na);
  /* bonds including H */
  if (h==AMBER_BOND_ALL || h==AMBER_BOND_WITH_H)
    for (i=0; i<nh; i++) {
      id1 = abs(t->bond_h[3*i]/3);
      id2 = abs(t->bond_h[3*i+1]/3);
      if (id1>=a0 && id2>=a0 && id1<=a1 && id2<=a1) {
        b[n_bonds][0] = id1;
        b[n_bonds][1] = id2;
        d[n_bonds] = t->bond_eq_val[t->bond_h[3*i+2]-1];
        n_bonds++;
        }
      }
  /* bonds without H */
  if (h==AMBER_BOND_ALL || h==AMBER_BOND_WITHOUT_H)
    for (i=0; i<na; i++) {
      id1 = abs(t->bond_ah[3*i]/3);
      id2 = abs(t->bond_ah[3*i+1]/3);
      if (id1>=a0 && id2>=a0 && id1<=a1 && id2<=a1) {
        b[n_bonds][0] = id1;
        b[n_bonds][1] = id2;
        d[n_bonds] = t->bond_eq_val[t->bond_ah[3*i+2]-1];
        n_bonds++;
        }
      }
  /* final size correction */
  if (n_bonds<(na+nh)) {
    b = mat_uresize_r(b,2,na+nh,n_bonds);
    d = vec_fresize(d,na+nh,n_bonds);
    }
  if (q)
    (*q) = d;
  else
    vec_ffree(d);
  (*n) = n_bonds;
  return(b);
  }

/* -------------------------------------------------------------------------- */
