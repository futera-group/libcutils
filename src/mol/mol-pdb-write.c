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

/* write one-atom line specification to file in pdb format
 
   atom_id   - atom sequential number (ID)
   atom_name - name of the atom (PDB format)
   res_id    - residuum sequential number (ID)
   res_name  - residuum name
   crd       - coordinates
   chrg      - charge
   prnt_chrg - print charge or not
   file      - output file stream */
void mol_pdb_fwrite_atom(unsigned atom_id, char *atom_name,
  unsigned res_id, char *res_name, double *crd, double chrg,
  short prnt_chrg, FILE *file) {
  /*  1 -  6 (6) ... keyword (record name) */
  fprintf(file,"%-6.6s","ATOM");
  /*  7 - 11 (5) ... serial number (atom number) */
  if (atom_id>99999)
    fprintf(file,"%05d",atom_id%100000);
  else
    fprintf(file,"%5d",atom_id);
  /* 12 - 12 (1) ... free space */
  fprintf(file," ");
  /* 13 - 16 (4) ... atom name */
  fprintf(file,"%-4.4s",atom_name);
  /* 17 - 17 (1) ... alternate location indicator */
  fprintf(file," ");
  /* 18 - 20 (3) ... residuum name */
  fprintf(file,"%-3.3s",res_name);
  /* 21 - 21 (1) ... free space */
  fprintf(file," ");
  /* 22 - 22 (1) ... chain identifier */
  fprintf(file," ");
  /* 23 - 26 (4) ... residuum sequence number */
  if (res_id>9999)
    fprintf(file,"%04d",res_id%10000);
  else
    fprintf(file,"%4d",res_id);
  /* 27 - 27 (1) ... residuum insertion code */
  fprintf(file," ");
  /* 28 - 30 (3) ... free space */
  fprintf(file,"   ");
  /* 31 - 38 (8) ... x-coordinate */
  fprintf(file,"%8.3f",crd[0]);
  /* 39 - 46 (8) ... y-coordinate */
  fprintf(file,"%8.3f",crd[1]);
  /* 47 - 54 (8) ... z-coordinate */
  fprintf(file,"%8.3f",crd[2]);
  /* 55 - 60 (6) ... occupancy */
  /* 61 - 66 (6) ... temperature factor */
  /* 77 - 78 (2) ... right-justified element symbol */
  /* 79 - 80 (2) ... charge */
  if (prnt_chrg)
    fprintf(file,"%8.4f",chrg);
  fprintf(file,"\n");
  }

/* write one-residuum data to file in pdb format

   r    - the residuum data
   file - output file stream */
void mol_pdb_fwrite_res(struct pdb_res *r, FILE *file) {
  unsigned i;
  for (i=0; i<r->n_atoms; i++)
    mol_pdb_fwrite_atom(r->atom[i].id,r->atom[i].name,r->id,r->name,
      r->atom[i].coord,r->atom[i].charge,1,file);
  }

/* write molecular data to file in pdb format
 
   p    - pointer to pdb molecular data struct 
   file - pointer to open file */
void mol_pdb_fwrite(struct pdb_mol *p, FILE *file) {
  unsigned i;
  if (!p)
    return;
  if (p->title && p->title[0]!='\0')
    fprintf(file,"TITLE  %-.73s\n",p->title);
  if (p->remark && p->remark[0]!='\0')
    fprintf(file,"REMARK %-.73s\n",p->remark);
  for (i=0; i<p->n_res; i++) {
    mol_pdb_fwrite_res(p->res+i,file);
    if (p->res[i].ter)
      fprintf(file,"TER\n");
    }
  fprintf(file,"END\n");
  }

/* write molecular data to file in pdb format
 
   p    - pointer to pdb molecular data struct 
   name - name of the file */
void mol_pdb_write(struct pdb_mol *p, char *name) {
  FILE *file = stdout;
  if (name && name[0])
    file = file_open(name,"w");
  mol_pdb_fwrite(p,file);
  if (name && name[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
