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
#include "prg/gauss.h"

/* -------------------------------------------------------------------------- */

/* merge two molecular orbitals data structs
 
   m1,m2 - the molecular orbital data structs
   n1,n2 - number of the molecular orbitals in each struct
   r1,r2 - number of the molecular orbitals coefficients in each struct */
struct gauss_mo* gauss_mo_merge(struct gauss_mo *m1, unsigned n1, unsigned r1,
  struct gauss_mo *m2, unsigned n2, unsigned r2) {
  unsigned i,j;
  struct gauss_mo *m;
  m = gauss_mo_new(n1+n2);
  for (i=0; i<n1; i++) {
    m[i].energy = m1[i].energy;
    m[i].coeff = vec_falloc(r1+r2);
    for (j=0; j<r1; j++)
      m[i].coeff[j] = m1[i].coeff[j];
    for (j=r1; j<(r1+r2); j++)
      m[i].coeff[j] = 0.0;
    }
  for (i=0; i<n2; i++) {
    m[n1+i].energy = m2[i].energy;
    m[n1+i].coeff = vec_falloc(r1+r2);
    for (j=0; j<r1; j++)
      m[n1+i].coeff[j] = 0.0;
    for (j=0; j<r2; j++)
      m[n1+i].coeff[r1+j] = m2[i].coeff[j];
    }
  return(m);
  }

/* -------------------------------------------------------------------------- */
