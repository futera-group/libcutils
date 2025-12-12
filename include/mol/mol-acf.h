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

#ifndef ZF_LIB_MOL_ACF_H
#define ZF_LIB_MOL_ACF_H

#include <stdio.h>

/* ACF atom struct */
struct acf_atom {
  char *name;               /* name of the atom */
  char *type;               /* type of the atom */
  unsigned num;             /* atomic number */
  unsigned res;             /* residuum number */
  double coord[3];          /* cartesian coordinates */
  double charge;            /* atomic charge */
  };

/* ACF molecule struct */
struct acf_mol {
  char *form;               /* chemical formula */
  char *name;               /* name of residuum */
  double charge;            /* molecular net charge */
  unsigned n_atoms;         /* number of atoms */
  unsigned n_bonds;         /* number of bonds */
  unsigned **bond;          /* bond array */
  struct acf_atom *atom;    /* atom array */
  };

/* allocate new acf molecule struct */
struct acf_mol *mol_acf_new(void);
/* allocate new acf atom data */
struct acf_atom *mol_acf_atom_new(unsigned);

/* free memory allocated for acf molecule struct */
void mol_acf_free(struct acf_mol*);
/* free memory allocated for acf atom data */
void mol_acf_atom_free(struct acf_atom*, unsigned);

/* coppy acf atom data */
void mol_acf_atom_copy(struct acf_atom*, struct acf_atom*);

/* read molecular data from acf file */
void mol_acf_read(struct acf_mol*, char*);
/* read molecular data from open acf file */
void mol_acf_fread(struct acf_mol*, FILE*);

/* write molecular data to file in acf format */
void mol_acf_write(struct acf_mol*, char*);
/* write molecular data to file in acf format */
void mol_acf_fwrite(struct acf_mol*, FILE*);

/* set empirical formula of molecular system in acf format */
void mol_acf_set_formula(struct acf_mol*);
/* calculate bonds for molecule in acf format */
void mol_acf_set_bonds(struct acf_mol*);

/* Try to assign Amber GAFF types to all atoms */
void mol_acf_type_gaff(struct acf_mol*);

/* check fragmentation and define additional bonds to connect all atoms */
void mol_acf_connect(struct acf_mol*);

/* convert acf molecular data struct to apc struct */
struct apc_mol *mol_acf_apc(struct acf_mol*);
/* convert acf molecular data struct to cif struct */
struct cif_mol *mol_acf_cif(struct acf_mol*, short, double*, double*);
/* convert acf molecular data struct to pdb struct */
struct pdb_mol *mol_acf_pdb(struct acf_mol*);
/* convert acf molecular data struct to xyz struct */
struct xyz_mol *mol_acf_xyz(struct acf_mol*);

#endif
