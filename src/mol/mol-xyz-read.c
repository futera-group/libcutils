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
#include <cmn/message.h>
#include <cmn/string.h>
#include "mol/atom.h"
#include "mol/molec.h"

/* -------------------------------------------------------------------------- */

/* read molecular data from xyz file
 
   x    - pointer to xyz molecular data struct 
   code - specify if terminate program in case of error or not
   file - pointer to open of the file */
int mol_xyz_fread(struct xyz_mol *x, short code, FILE *file) {
  unsigned i;
  char sym[80];
  char *line;
  /* number of atoms */
  line = str_read_line_new(file);
  if (!line) {
    if (code)
      msg_error("xyz file is empty",code);
    else
      return(0);
    }
  if (sscanf(line,"%u",&(x->n_atoms))!=1) {
    if (code)
      msg_error("cannot read number of atoms in xyz file",code);
    else
      return(0);
    }
  str_free(line);
  /* title */
  x->title = str_read_line_new(file);
  if (!x->title) {
    if (code)
      msg_error("cannot read title line in xyz file",code);
    else
      return(0);
    }
  str_trim(x->title);
  /* re/allocate memory for atom data */
  mol_xyz_atom_free(x->atom,x->n_atoms);
  x->atom = mol_xyz_atom_new(x->n_atoms);
  /* read atom coordinates */
  for (i=0; i<x->n_atoms; i++) {
    line = str_read_line_new(file);
    if (!line) {
      if (code)
        msg_error("unexpected end of xyz file while reading atoms",code);
      else
        return(0);
      }
    if (sscanf(line,"%s%lf%lf%lf",sym,&(x->atom[i].coord[0]),
      &(x->atom[i].coord[1]),&(x->atom[i].coord[2]))!=4) {
      if (code)
        msg_error("invalid format of xyz file",code);
      else
        return(0);
      }
    x->atom[i].name = str_copy_new(sym);
    x->atom[i].num = atom_num(sym);
    str_free(line);
    }
  return(1);
  }

/* read molecular data from xyz file
 
   x    - pointer to xyz molecular data struct 
   name - name of the file */
void mol_xyz_read(struct xyz_mol *x, char *name) {
  FILE *file = file_open(name,"r");
  mol_xyz_fread(x,1,file);
  file_close(file);
  }

/* -------------------------------------------------------------------------- */
