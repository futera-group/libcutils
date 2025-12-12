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

#include <cmn/vector.h>
#include "qmc/basis.h"

/* -------------------------------------------------------------------------- */

/* copy data from one basis set shell data struct to another
 
   c1 - the destination data struct where the data are writen
   c2 - the source data struct where the data are read */
void basis_shell_copy(struct basis_shell *s1, struct basis_shell *s2) {
  s1->type = s2->type;
  s1->bf = s2->bf;
  s1->n_bfce = s2->n_bfce;
  s1->n_prim = s2->n_prim;
  s1->exp = vec_fcopy_new(s2->exp,s2->n_prim);
  s1->cf1 = vec_fcopy_new(s2->cf1,s2->n_prim);
  s1->cf2 = vec_fcopy_new(s2->cf2,s2->n_prim);
  }

/* create copy of array of basis set shell data structure
 
   s0 - the source data struct where the data are read
   n  - number of the shells */
struct basis_shell* basis_shell_copy_new(struct basis_shell *s0, unsigned n) {
  unsigned i;
  struct basis_shell *s1;
  s1 = basis_shell_new(n);
  for (i=0; i<n; i++) {
    s1[i].type = s0[i].type;
    s1[i].bf = s0[i].bf;
    s1[i].n_bfce = s0[i].n_bfce;
    s1[i].n_prim = s0[i].n_prim;
    s1[i].exp = vec_fcopy_new(s0[i].exp,s0[i].n_prim);
    s1[i].cf1 = vec_fcopy_new(s0[i].cf1,s0[i].n_prim);
    s1[i].cf2 = vec_fcopy_new(s0[i].cf2,s0[i].n_prim);
    }
  return(s1);
  }

/* copy data from one basis set center data struct to another
 
   c1 - the destination data struct where the data are writen
   c2 - the source data struct where the data are read */
void basis_center_copy(struct basis_center *c1, struct basis_center *c2) {
  unsigned i;
  c1->bf = c2->bf;
  c1->type = c2->type;
  c1->n_bfce = c2->n_bfce;
  c1->n_shells = c2->n_shells;
  c1->shell = basis_shell_new(c2->n_shells);
  for (i=0; i<3; i++)
    c1->coord[i] = c2->coord[i];
  for (i=0; i<c1->n_shells; i++)
    basis_shell_copy(c1->shell+i,c2->shell+i);
  }

/* copy array of basis set center data structures
 
   c0 - the source data struct where the data are read
   n  - number of the centers */
struct basis_center* basis_center_copy_new(struct basis_center *c0,
  unsigned n) {
  unsigned i,j;
  struct basis_center *c1;
  c1 = basis_center_new(n);
  for (i=0; i<n; i++) {
    c1[i].bf = c0[i].bf;
    c1[i].type = c0[i].type;
    c1[i].n_bfce = c0[i].n_bfce;
    c1[i].n_shells = c0[i].n_shells;
    c1[i].shell = basis_shell_copy_new(c0[i].shell,c0[i].n_shells);
    for (j=0; j<3; j++)
      c1[i].coord[j] = c0[i].coord[j];
    }
  return(c1);
  }

/* -------------------------------------------------------------------------- */
