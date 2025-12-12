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
#include <cmn/message.h>
#include <cmn/string.h>
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* print list of residui or atoms

   p    - pointer to pdb molecular data struct
   code - specify what to print
   file - pointer to open file where to print */
void mol_pdb_print(struct pdb_mol *p, short code, FILE *file) {
  unsigned i,j,n = 0,row = 1;
  if (!p)
    return;
  /* print header */
  switch (code) {
    case PDB_PRINT_RES:
    case PDB_PRINT_RES_ALL:
      fprintf(file,"\n  Residues:\n\n   "); break;
    case PDB_PRINT_ATOM:
    case PDB_PRINT_ATOM_ALL: 
      fprintf(file,"\n  Atoms:\n\n   "); break;
    default: msg_error("invalid specification for PDB data printing",1);
    }
  for (i=0; i<10; i++)
    fprintf(file,"%5d",(i+1)%10); /* print data */
  for (i=0; i<p->n_res; i++) 
    for (j=0; j<p->res[i].n_atoms; j++) {
      /* skip water if required */
      if ((code==PDB_PRINT_RES || code==PDB_PRINT_ATOM) &&
        (str_compare(p->res[i].name,"WAT") ||
         str_compare(p->res[i].name,"HOH")))
        continue;
      /* print info */
      else {
        /* residuum name */
        if (code==PDB_PRINT_RES || code==PDB_PRINT_RES_ALL) {
          if (n%10==0) {
            fprintf(file,"\n%5d",row);
            row += 10;
            }
          fprintf(file,"%5s",p->res[i].name);
          n++;
          break;
          }
        /* atom name */
        else {
          if (n%10==0) {
            fprintf(file,"\n%5d",row);
            row += 10;
            }
          fprintf(file,"%5s",p->res[i].atom[j].name);
          n++;
          }
        }
      }
  fprintf(file,"\n");
  }

/* -------------------------------------------------------------------------- */
