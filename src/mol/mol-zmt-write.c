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
#include <cmn/queue.h>
#include <cmn/string.h>
#include "mol/atom.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* write molecular data to file in zmt format 
 
   z    - pointer to zmt molecular data struct 
   frmt - format of z-matrix 
   file - pointer to open file */
void mol_zmt_fwrite(struct zmt_mol *z, short frmt, FILE *file) {
  char sym[256];
  unsigned i,n,nb,na,nd;
  struct queue *q = NULL;
  nb = na = nd = 1;
  n = str_unum_len(z->n_atoms);
  if (frmt==ZMT_FRMT_VAR)
    q = queue_alloc();
  /* print out z-matrix */
  for (i=0; i<z->n_atoms; i++) {
    fprintf(file," %-3s",atom_name(z->atom[i].num));
    /* bonds */
    if (i>0) {
      fprintf(file," %3d",z->atom[i].bond_id+1);
      if (frmt==ZMT_FRMT_VAR) {
        sprintf(sym,"b%0*d",n,nb++);
        fprintf(file," %-6s",sym);
        queue_sadd(q,sym);
        }
      else
        fprintf(file," %8.3f",z->atom[i].bond_val);
      }
    /* angles */
    if (i>1) {
      fprintf(file," %3d",z->atom[i].angle_id+1);
      if (frmt==ZMT_FRMT_VAR) {
        sprintf(sym,"a%0*d",n,na++);
        fprintf(file," %-6s",sym);
        queue_sadd(q,sym);
        }
      else
        fprintf(file," %8.3f",z->atom[i].angle_val);
      }
    /* dihedral angles */
    if (i>2) {
      fprintf(file," %3d",z->atom[i].dihed_id+1);
      if (frmt==ZMT_FRMT_VAR) {
        sprintf(sym,"d%0*d",n,nd++);
        fprintf(file," %-6s",sym);
        queue_sadd(q,sym);
        }
      else 
        fprintf(file," %8.3f",z->atom[i].dihed_val);
      }
    fprintf(file,"\n");
    }
  if (frmt==ZMT_FRMT_VAR) {
    /* print out internal coordinate values */
    fprintf(file,"\n");
    for (i=0; i<z->n_atoms; i++) {
      if (i>0) {
        queue_sget(q,sym);
        fprintf(file,"%-12s %12.5f\n",sym,z->atom[i].bond_val);
        }
      if (i>1) {
        queue_sget(q,sym);
        fprintf(file,"%-12s %12.5f\n",sym,z->atom[i].angle_val);
        }
      if (i>2) {
        queue_sget(q,sym);
        fprintf(file,"%-12s %12.5f\n",sym,z->atom[i].dihed_val);
        }
      }
    /* clean memory */
    queue_free(q);
    }
  }

/* write molecular data to file in zmt format
 
   x    - pointer to zmt molecular data struct 
   fmrt - format of z-matrix 
   name - name of the file */
void mol_zmt_write(struct zmt_mol *x, short frmt, char *name) {
  FILE *file = stdout;
  if (name && name[0])
    file = file_open(name,"w");
  mol_zmt_fwrite(x,frmt,file);
  if (name && name[0])
    file_close(file);
  }

/* -------------------------------------------------------------------------- */
