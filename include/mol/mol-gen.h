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

#ifndef ZF_LIB_MOL_GEN_H
#define ZF_LIB_MOL_GEN_H

#include "mol/cell.h"

/* general atom struct */
struct gen_atom {
  unsigned num;                        /* atomic number */
  double coord[3];                     /* atomic coordinates */
  double charge;                       /* charge */
  };

/* general molecule struct */
struct gen_mol {
  char *title;                         /* title line */
  unsigned n_atoms;                    /* number of atoms */
  struct cell *cell;                   /* unit cell */
  struct gen_atom *atom;               /* array of atoms */
  };

/* allocate new general molecule struct */
struct gen_mol *mol_gen_new(void);
/* allocate array of general atom data structs */
struct gen_atom* mol_gen_atom_new(unsigned);

/* free memory allocated for general molecule struct */
void mol_gen_free(struct gen_mol*);
/* free memory allocated for general atom data */
void mol_gen_atom_free(struct gen_atom*);

/* create bond array with atom indices for the given molecule */
unsigned **mol_gen_bonds(struct gen_mol*, unsigned*);

/* convert generic molecular data struct to cif format */
struct cif_mol* mol_gen_cif(struct gen_mol*);
/* convert generic molecular data struct to xyz format */
struct xyz_mol* mol_gen_xyz(struct gen_mol*);

#endif
