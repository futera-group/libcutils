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

#include "prg/amber.h"

/* -------------------------------------------------------------------------- */

/* wrap atom coordinates into the PBC box, keep compact residui
 
   t - pointer to amber topology struct
   x - array with coordinates
   b - PBC box data */
void amber_crd_wrap(struct amber_top *t, double *x, double *b) {
  unsigned i,j,r,a0,a1;
  double c[3];
  int s;
  for (r=0; r<t->pointers[AMBER_POINTER_NRES]; r++) {
    amber_res_mass_center(t,x,r,c);
    amber_res_atom_id(t,r,&a0,&a1);
    for (i=0; i<3; i++)
      if (c[i]<0.0 || c[i]>b[i]) {
        s = (int)(c[i]/b[i]);
        if (c[i]<0.0)
          s--;
        for (j=a0; j<a1; j++)
          x[3*j+i] -= (s*b[i]);
        }
    }
  }

/* -------------------------------------------------------------------------- */
