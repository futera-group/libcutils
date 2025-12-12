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

#ifndef ZF_LIB_MOL_APC_H
#define ZF_LIB_MOL_APC_H

#include <stdio.h>

/* APC atom struct */
struct apc_atom {
  char *name;                 /* name of the atom */
  char *type;                 /* type of the atom */
  unsigned tree;              /* tree mark */
  unsigned id;                /* atom ID */
  double coord[3];            /* atomic coordinates */
  double charge;              /* atomic charge */
  };

/* APC molecule struct */
struct apc_mol {
  char *title;                /* title line */
  char *file;                 /* file name */
  char *resname;              /* file name */
  unsigned n_atoms;           /* number of atoms */
  unsigned n_loops;           /* number of loops */
  unsigned n_imprs;           /* number of improper dihedrals */
  unsigned **loop;            /* loop array */
  unsigned **impr;            /* improper array */
  struct apc_atom *atom;      /* atom array */
  };

/* allocate new apc molecule struct */
struct apc_mol *mol_apc_new(void);
/* allocate array of apc atom data structures */
struct apc_atom *mol_apc_atom_new(unsigned);

/* free memory allocated for apc molecule struct */
void mol_apc_free(struct apc_mol*);
/* free memory allocater for array of apc atom data structures */
void mol_apc_atom_free(struct apc_atom*, unsigned);

/* copy apc atom data */
void mol_apc_atom_copy(struct apc_atom*, struct apc_atom*);

/* read molecular data from apc file */
void mol_apc_read(struct apc_mol*, char*);
/* read molecular data from open apc file */
void mol_apc_fread(struct apc_mol*, FILE*);

/* write molecular data to file in apc format */
void mol_apc_write(struct apc_mol*, short, short, char*);
/* write molecular data to file in apc format */
void mol_apc_fwrite(struct apc_mol*, short, short, FILE*);

/* return ID of the atom specified by name */
unsigned mol_apc_atom_id(struct apc_mol*, char*);
/* return ID of specified atom */
unsigned mol_apc_atom_id2id(struct apc_mol*, unsigned);
/* reset all atom IDs as well as loop and improper pointers */
void mol_apc_atom_id_set(struct apc_mol*);

/* convert tree mark to internal code */
short mol_apc_tree_id(char*);
/* convert tree internal code to tree mark */
void mol_apc_tree_mark(short, char*);

/* convert apc molecular data struct to acf struct */
struct acf_mol *mol_apc_acf(struct apc_mol*);
/* convert apc molecular data struct to cif struct */
struct cif_mol *mol_apc_cif(struct apc_mol*, short, double*, double*);
/* convert apc molecular data struct to pdb struct */
struct pdb_mol *mol_apc_pdb(struct apc_mol*);
/* convert apc molecular data struct to xyz struct */
struct xyz_mol *mol_apc_xyz(struct apc_mol*);

#endif
