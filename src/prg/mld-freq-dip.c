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

#include "cmn/vector.h"
#include "prg/molden.h"

/* -------------------------------------------------------------------------- */

/* estimate electric transition dipole moment from the normal mode displacement
 
   a - array of atoms
   d - array of atomic displacements
   n - number of atoms
   m - dipole coordinates (output) */
void mld_freq_dip(struct mld_atom *a, double **d, unsigned n, double *m) {
  double q,m0[3],m1[3];
  unsigned i,j;
  /* charge */
  q = 0.0;
  for (i=0; i<n; i++)
    q += a[i].num;
  q /= n;
  /* dipole moment of the relaxed structure */
  vec_fset(m0,0.0,3);
  for (i=0; i<n; i++)
    for (j=0; j<3; j++)
      m0[j] += ((a[i].num - q) * a[i].crd[j]);
  /* dipole moment of the displaced structure */
  vec_fset(m1,0.0,3);
  for (i=0; i<n; i++)
    for (j=0; j<3; j++)
      m1[j] += ((a[i].num - q) * (a[i].crd[j] + d[i][j]));
  /* transition dipole moment */
  for (i=0; i<3; i++)
    m[i] = m1[i] - m0[i];
  }

/* -------------------------------------------------------------------------- */
