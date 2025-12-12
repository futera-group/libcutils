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
#include <cmn/file.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* write molecular data to file in acf format
 
   c - pointer to acf molecular data struct
   f - pointer to open acf file */
void mol_acf_fwrite(struct acf_mol *c, FILE *f) {
  unsigned i;
  fprintf(f,"CHARGE%10.2f ( %d )\n",c->charge,(int)(c->charge));
  fprintf(f,"Formula: %s\n",c->form);
  for (i=0; i<c->n_atoms; i++) 
    fprintf(f,"ATOM %6d %-4s %3s%6d%12.3f%8.3f%8.3f%10.6f%10s\n",i+1,
      c->atom[i].name,c->name,1,c->atom[i].coord[0],c->atom[i].coord[1],
      c->atom[i].coord[2],c->atom[i].charge,c->atom[i].type);
  for (i=0; i<c->n_bonds; i++)
    fprintf(f,"BOND %4d %4d %4d %4d%8s%5s\n",i+1,
      c->bond[i][0]+1,c->bond[i][1]+1,1,c->atom[c->bond[i][0]].name,
      c->atom[c->bond[i][1]].name);
  }

/* write molecular data to file in acf format
 
   c - pointer to acf molecular data struct
   f - name of the file */
void mol_acf_write(struct acf_mol *c, char *f) {
  FILE *file = stdout;
  if (f && f[0])
    file = file_open(f,"w");
  mol_acf_fwrite(c,file);
  if (f && f[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
