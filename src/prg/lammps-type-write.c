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

#include <stdio.h>
#include "prg/lammps.h"

/* -------------------------------------------------------------------------- */

/* write atomic types to LAMMPS data file
 
   m - array of the types
   n - number of atomic types
   f - open file stream */
void lammps_type_write(struct lammps_type *m, unsigned n, FILE *f) {
  unsigned i;
  if (m && n) {
    fprintf(f,"\nMasses\n\n");
    for (i=0; i<n; i++) {
      fprintf(f," %5d %12.6f",i+1,m[i].mass);
      if (m[i].name[0])
        fprintf(f,"  # %s",m[i].name);
      fprintf(f,"\n");
      }
    }
  }

/* -------------------------------------------------------------------------- */
