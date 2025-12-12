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
#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

/* return ID of residuum for given atom
 
   t  - pointer to amber topology struct
   a  - ID of the atom in the residuum (output)
   id - system ID of the atom (input) */
unsigned amber_res_id(struct amber_top *t, unsigned *a, unsigned id) {
  unsigned i,a0,a1;
  for (i=0; i<t->pointers[AMBER_POINTER_NRES]; i++) {
    a0 = t->res_pointer[i];
    a1 = (i==(t->pointers[AMBER_POINTER_NRES]-1) ? 
      t->pointers[AMBER_POINTER_NATOM] : t->res_pointer[i+1]-1);
    if ((id+1)>=a0 && (id+1)<=a1) {
      if (a) 
        (*a) = id+1-a0;
      return(i);
      }
    }
  return(0);
  }

/* check if the given residuum is terminal of not
 
   t  - pointer to amber topology struct 
   id - ID of the residuum */
short amber_res_is_ter(struct amber_top *t, unsigned id) {
  unsigned i,j,af,al;
  af = t->res_pointer[id]-1;
  if (t->pointers[AMBER_POINTER_NRES]-1==id)
   al = t->pointers[AMBER_POINTER_NATOM]-1;
  else
   al = t->res_pointer[id+1]-2;
  for (i=af; i<=al; i++) {
    for (j=0; j<t->pointers[AMBER_POINTER_NBONH]; j++)
      if (t->bond_h[3*j]/3==i && t->bond_h[3*j+1]/3>al) {
        return(0);
        break;
        }
    for (j=0; j<t->pointers[AMBER_POINTER_NBONA]; j++)
      if (t->bond_ah[3*j]/3==i && t->bond_ah[3*j+1]/3>al) {
        return(0);
        break;
        }
    }
  return(1);
  }

/* return array with IDs of unique residui
 
   t - pointer to amber topology struct
   n - number of unique residui (output) */
unsigned *amber_res_unique(struct amber_top *t, unsigned *n) {
  unsigned i,j,nr=0,*r=NULL;
  short unique;
  r = vec_ualloc(t->pointers[AMBER_POINTER_NRES]);
  for (i=0; i<t->pointers[AMBER_POINTER_NRES]; i++) {
    unique = 1;
    for (j=0; j<nr; j++)
      if (str_compare(t->res_names[i],t->res_names[r[j]])) {
        unique = 0;
        break;
        }
    if (unique)
      r[nr++] = i;
    }
  r = vec_uresize(r,t->pointers[AMBER_POINTER_NRES],nr);
  (*n) = nr;
  return(r);
  }

/* return IDs of first and last+1 atom of specified residuum
 
   t  - pointer to amber topology struct
   a0 - first atom of the residuum
   a1 - last+1 atom of the residuum */
void amber_res_atom_id(struct amber_top *t, unsigned r,
  unsigned *a0, unsigned *a1) {
  if (a0)
    (*a0) = t->res_pointer[r]-1;
  if (a1)
    (*a1) = (r==(t->pointers[AMBER_POINTER_NRES]-1) ?
      t->pointers[AMBER_POINTER_NATOM] : t->res_pointer[r+1]-1);
  }

/* return number of atoms in specified residuum
 
   t - pointer to amber topology struct
   r - the residuum ID */
unsigned amber_res_natoms(struct amber_top *t, unsigned r) {
  unsigned a0,a1;
  amber_res_atom_id(t,r,&a0,&a1);
  return(a1-a0);
  }

/* calculate center of mass of specified residuum of the system

   t - pointer to amber topology struct
   x - array with coordinates
   r - the residuum ID
   c - array for COM coordinates */
void amber_res_mass_center(struct amber_top *t, double *x, unsigned r,
  double *c) {
  unsigned i,j,a0,a1;
  double m;
  amber_res_atom_id(t,r,&a0,&a1);
  m = c[0] = c[1] = c[2] = 0.0;
  for (i=a0; i<a1; i++) {
    for (j=0; j<3; j++)
      c[j] += (t->mass[i]*x[3*i+j]);
    m += t->mass[i];
    }
  for (i=0; i<3; i++)
    c[i] /= m;
  }

/* -------------------------------------------------------------------------- */
