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

#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* cyclic shift of data array

   c - pointer to xyz molecular data struct 
   s - direction indicator */
void mol_cub_data_shift(struct cub_mol *c, short s) {
  unsigned i,j,k;
  double v;
  switch (s) {
    case 0: 
      for (i=0; i<c->grid->range[1]; i++)
        for (j=0; j<c->grid->range[2]; j++) {
          v = c->data[c->grid->range[0]-1][i][j];
          for (k=c->grid->range[0]-1; k>0; k--)
            c->data[k][i][j] = c->data[k-1][i][j];
          c->data[0][i][j] = v;
          }
      break;
    case 1:
      for (i=0; i<c->grid->range[0]; i++)
        for (j=0; j<c->grid->range[2]; j++) {
          v = c->data[i][c->grid->range[1]-1][j];
          for (k=c->grid->range[1]-1; k>0; k--)
            c->data[i][k][j] = c->data[i][k-1][j];
          c->data[i][k][j] = v;
          }
      break;
    case 2:
      for (i=0; i<c->grid->range[0]; i++)
        for (j=0; j<c->grid->range[1]; j++) {
          v = c->data[i][j][c->grid->range[2]-1];
          for (k=c->grid->range[2]-1; k>0; k--)
            c->data[i][j][k] = c->data[i][j][k-1];
          c->data[i][j][k] = v;
          }
      break;
    }
  }

/* shift data in cube and apply periodic boundary conditions

   c - pointer to xyz molecular data struct
   s - shifting vector  */
void mol_cub_shift(struct cub_mol *c, double *s) {
  unsigned i,j,n[3];
  double dm[3];
  /* initialization */
  for (i=0; i<3; i++) {
    dm[i] = c->grid->range[i]*c->grid->vector[i][i];
    while (s[i]<0.0)
      s[i] += dm[i];
    while (s[i]>dm[i])
      s[i] -= dm[i];
    n[i] = s[i]/c->grid->vector[i][i];
    }
  /* shift structure coordinates */
  for (i=0; i<c->n_atoms;  i++)
    for (j=0; j<3; j++) {
      c->atom[i].coord[j] += s[j];
      if (c->atom[i].coord[j]>dm[j])
        c->atom[i].coord[j] -= dm[j];
      }
  /* shift grid data */
  for (i=0; i<3; i++)
    for (j=0; j<n[i]; j++)
      mol_cub_data_shift(c,i);
  }

/* -------------------------------------------------------------------------- */
