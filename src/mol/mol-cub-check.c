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

#include <math.h>
#include <stdio.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* check if two cube data structs are compatible 

   c1,c2 - pointer to cube data structs
   c     - error code in case of incompatibility
   m     - incompatibility message
   at    - check atom number and types */
short mol_cub_check_compatibility(struct cub_mol *c1, struct cub_mol *c2,
  short *c, char *m, short at) {
  unsigned i,j;
  /* atomic data */
  if (at) {
    if (c1->n_atoms!=c2->n_atoms) {
      if (m)
        sprintf(m,"different number of atoms: %d vs %d",c1->n_atoms,c2->n_atoms);
      if (c)
        (*c) = CUB_ICMP_NATOM;
      return(0);
      }
    for (i=0; i<c1->n_atoms; i++)
      if (c1->atom[i].num!=c2->atom[i].num) {
        if (m)
          sprintf(m,"different atomic types of atom %d: %d vs %d",
            i+1,c1->atom[i].num,c2->atom[i].num);
        if (c)
          (*c) = CUB_ICMP_ATYPE;
        return(0);
        }
    }
  /* grid size */
  for (i=0; i<3; i++)
    if (c1->grid->range[i]!=c2->grid->range[i]) {
      if (m)
        sprintf(m,"different grid size: %d-%d-%d vs %d-%d-%d",
          c1->grid->range[0],c1->grid->range[1],c1->grid->range[2],
          c2->grid->range[0],c2->grid->range[1],c2->grid->range[2]);
      if (c)
        (*c) = CUB_ICMP_GRID;
      return(0);
      }
  /* lattice specification */
  for (i=0; i<3; i++)
    if (fabs(c1->grid->origin[i]-c2->grid->origin[i])>CUB_CHECK_THRD) {
      if (m)
        sprintf(m,"different lattice origin:"
          " [%.6f,%.6f,%.6f] vs [%.6f,%.6f,%.6f]",
          c1->grid->origin[0],c1->grid->origin[1],c1->grid->origin[2],
          c2->grid->origin[0],c2->grid->origin[1],c2->grid->origin[2]);
      if (c)
        (*c) = CUB_ICMP_LORIG;
      return(0);
      }
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      if (fabs(c1->grid->vector[i][j]-c2->grid->vector[i][j])>CUB_CHECK_THRD) {
        if (m)
        sprintf(m,"different lattice vectors %d:"
          " [%.6f,%.6f,%.6f] vs [%.6f,%.6f,%.6f]",i+1,
          c1->grid->vector[i][0],c1->grid->vector[i][1],c1->grid->vector[i][2],
          c2->grid->vector[i][0],c2->grid->vector[i][1],c2->grid->vector[i][2]);
        if (c)
          (*c) = CUB_ICMP_LVEC;
        return(0);
        }
  return(1);
  }

/* -------------------------------------------------------------------------- */
