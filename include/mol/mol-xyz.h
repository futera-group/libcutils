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

#ifndef ZF_LIB_MOL_XYZ_H
#define ZF_LIB_MOL_XYZ_H

#include <stdio.h>

/* XYZ atom struct */
struct xyz_atom {
  char *name;               /* atomic name (PDB) */
  unsigned num;             /* atomic number */
  double coord[3];          /* atomic coordinates */
  };

/* XYZ molecule struct */
struct xyz_mol {
  char *title;              /* title line */
  unsigned n_atoms;         /* number of atoms */
  struct xyz_atom *atom;    /* vector of atoms */
  };

/* allocate new xyz molecule struct */
struct xyz_mol *mol_xyz_new(void);
/* allocate array of xyz atom data struct */
struct xyz_atom *mol_xyz_atom_new(unsigned);

/* free memory allocated for xyz molecule struct */
void mol_xyz_free(struct xyz_mol*);
/* free memory allocated for xyz atom array */
void mol_xyz_atom_free(struct xyz_atom*, unsigned);

/* read molecular data from xyz file */
void mol_xyz_read(struct xyz_mol*, char*);
/* read molecular data from open xyz file */
int mol_xyz_fread(struct xyz_mol*, short, FILE*);

/* write molecular data to file in xyz format */
void mol_xyz_write(struct xyz_mol*, char*);
/* write molecular data to file in xyz format */
void mol_xyz_fwrite(struct xyz_mol*, FILE*);
/* write molecular data to file in xyz format with specific number format */
void mol_xyz_fwritef(struct xyz_mol*, char*, FILE*);

/* copy XYZ molecular data struct */
struct xyz_mol *mol_xyz_copy_new(struct xyz_mol*);

/* merge two XYZ molecular data structs */
struct xyz_mol *mol_xyz_merge(struct xyz_mol*, struct xyz_mol*);

/* create bond array with atom indices for the given molecule */
unsigned **mol_xyz_bonds(struct xyz_mol*, unsigned*);

/* delete specified atom from the structure */
void mol_xyz_atom_del(struct xyz_mol*, unsigned);
/* Delete atom with specified name from the structure */
void mol_xyz_atom_del_name(struct xyz_mol*, char*);
/* order atoms according to new indices */
void mol_xyz_atom_order(struct xyz_mol*, unsigned*);
/* find atom in the xyz data structure and return its id */
unsigned mol_xyz_atom_find(struct xyz_mol*, char*, short*);

/* set coordinates to XYZ data struct from general molecular struct */
void mol_xyz_set_coord_gen(struct xyz_mol*, struct gen_mol*);

/* convert xyz molecular data struct to acf struct */
struct acf_mol *mol_xyz_acf(struct xyz_mol*);
/* convert xyz molecular data struct to apc struct */
struct apc_mol *mol_xyz_apc(struct xyz_mol*);
/* convert xyz molecular data struct to cif struct */
struct cif_mol *mol_xyz_cif(struct xyz_mol*, short, double*, double*);
/* convert xyz molecular data struct to gen struct */
struct gen_mol *mol_xyz_gen(struct xyz_mol*);
/* convert xyz molecular data struct to pdb struct */
struct pdb_mol *mol_xyz_pdb(struct xyz_mol*, short);
/* convert xyz struct to zmt molecular data struct */
struct zmt_mol *mol_xyz_zmt(struct xyz_mol*, struct zmt_mol*);

#endif
