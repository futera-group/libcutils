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

#ifndef ZF_LIB_MOL_CIF_H
#define ZF_LIB_MOL_CIF_H

#include <stdio.h>
#include "mol/cell.h"

/* CIF atom struct */
struct cif_atom {
  unsigned num;                /* atomic number */
  char *name;                  /* atom label */
  double coord[3];             /* atomic coordinates */
  };

/* CIF molecule struct */
struct cif_mol {
  unsigned n_atoms;            /* number of atoms */
  struct cell *cell;           /* simulation cell */
  struct cif_atom *atom;       /* vector of atoms */
  };

/* allocate new cif molecule struct */
struct cif_mol* mol_cif_new(void);
/* allocate array of cif atom data struct */
struct cif_atom* mol_cif_atom_new(unsigned);

/* free memeory allocated for cif molecule struct */
void mol_cif_free(struct cif_mol*);
/* free memory allocated for array of cif atom data structures */
void mol_cif_atom_free(struct cif_atom*, unsigned);

/* copy cif molecular data struct */
struct cif_mol* mol_cif_copy_new(struct cif_mol*);
/* copy cif atom data */
void mol_cif_atom_copy(struct cif_atom*, struct cif_atom*);

/* read molecular data from cif file */
void mol_cif_read(struct cif_mol*, char*);
/* read molecular data from already open cif file */
void mol_cif_fread(struct cif_mol*, FILE*);

/* write molecular data to cif file */
void mol_cif_write(struct cif_mol*, char*);
/* write molecular data to already open cif file */
void mol_cif_fwrite(struct cif_mol*, FILE*);

/* create bond array with atom indices for the given molecule */
unsigned **mol_cif_bonds(struct cif_mol*, unsigned*);

/* convert cif molecular data struct to acf struct */
struct acf_mol *mol_cif_acf(struct cif_mol*);
/* convert cif molecular data struct to apc struct */
struct apc_mol *mol_cif_apc(struct cif_mol*);
/* convert cif molecular data struct to pdb struct */
struct pdb_mol *mol_cif_pdb(struct cif_mol*);
/* convert cif molecular data struct to xyz struct */
struct xyz_mol *mol_cif_xyz(struct cif_mol*);

#endif
