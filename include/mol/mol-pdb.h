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

#ifndef ZF_LIB_MOL_PDB_H
#define ZF_LIB_MOL_PDB_H

#include <stdio.h>

#define PDB_PRINT_RES       1  /* print all residui names but water */
#define PDB_PRINT_RES_ALL   2  /* print all residui names */
#define PDB_PRINT_ATOM      3  /* print all atom names but water */
#define PDB_PRINT_ATOM_ALL  4  /* print all atom names */

/* PDB atom struct */
struct pdb_atom {
  char *name;               /* name of the atom */
  unsigned id;              /* atom ID */
  double coord[3];          /* cartesian coordinates */
  double charge;            /* atomic charge */
  };

/* PDB residuum struct */
struct pdb_res {
  char *name;               /* name of the residuum */
  short ter;                /* TER flag */
  unsigned id;              /* residuum ID */
  unsigned n_atoms;         /* number of atoms */
  struct pdb_atom *atom;    /* vector of atoms */
  };

/* PDB molecule struct */
struct pdb_mol {
  char *title;              /* title line */
  char *remark;             /* remark line */
  unsigned n_res;           /* number of residui */
  struct pdb_res *res;      /* vector of residui */
  };

/* PDB atom struct - first read */
struct pdb_unit {
  short ter;                /* terminating atom */
  char *name;               /* name of the atom */
  char *resname;            /* name of the residuum */
  unsigned resid;           /* residuum number */
  double coord[3];          /* atomic coordinates */
  double charge;            /* atomic charge */
  };

/* allocate new pdb molecule struct */
struct pdb_mol *mol_pdb_new(void);
/* allocate memory for array of residui */
struct pdb_res *mol_pdb_res_new(unsigned);
/* allocate memory for array of atoms */
struct pdb_atom *mol_pdb_atom_new(unsigned);
/* allocate memory for new temporary info pdb atom struct */
struct pdb_unit *mol_pdb_unit_new(void);

/* free memory allocated for pdb molecule struct */
void mol_pdb_free(struct pdb_mol*);
/* free memory allocated for array of pdb residui */
void mol_pdb_res_free(struct pdb_res*, unsigned);
/* free memory allocated for array of pdb atoms */
void mol_pdb_atom_free(struct pdb_atom*, unsigned);
/* free memory allocated for pbc unit data structure */
void mol_pdb_unit_free(struct pdb_unit*);

/* read molecular data from pdb file */
void mol_pdb_read(struct pdb_mol*, char*);
/* read molecular data from open pdb file */
int mol_pdb_fread(struct pdb_mol*, short, FILE*);

/* write molecular data to file in pdb format */
void mol_pdb_write(struct pdb_mol*, char*);
/* write molecular data to file in pdb format */
void mol_pdb_fwrite(struct pdb_mol*, FILE*);
/* write one-residuum data to file in pdb format */
void mol_pdb_fwrite_res(struct pdb_res*, FILE*);
/* write one-atom line specification to file in pdb format */
void mol_pdb_fwrite_atom(unsigned, char*, 
  unsigned, char*, double*, double, short, FILE*);

/* return number of atoms in pdb struct */
unsigned mol_pdb_atom_n(struct pdb_mol*);
/* find residuum and atom ID for given atom */
short mol_pdb_atom_id(struct pdb_mol*, unsigned, unsigned*, unsigned*);
/* set unique pdb atomic name using the sequential number */
void mol_pdb_atom_name_id(char*, unsigned, unsigned);
/* set pdb atomic name and fix the letter position */
void mol_pdb_atom_name_pos(struct pdb_mol*, unsigned, unsigned, char*);

/* find atom in specified residuum and returns its ID */
unsigned mol_pdb_res_atom_find(struct pdb_res*, char*, short*);

/* create bond array with atom indices for the given molecule */
unsigned **mol_pdb_bonds(struct pdb_mol*, unsigned, unsigned*);

/* merge two PDB molecular data structs */
struct pdb_mol *mol_pdb_merge(struct pdb_mol*, struct pdb_mol*);

/* print list of residui or atoms */
void mol_pdb_print(struct pdb_mol*, short, FILE*);

/* set coordinates to PDB data struct from general molecular struct */
void mol_pdb_set_coord_gen(struct pdb_mol*, struct gen_mol*);

/* convert pdb molecular data struct to acf struct */
struct acf_mol *mol_pdb_acf(struct pdb_mol*);
/* convert pdb molecular data struct to apc struct */
struct apc_mol *mol_pdb_apc(struct pdb_mol*, unsigned);
/* convert pdb molecular data struct to cif struct */
struct cif_mol *mol_pdb_cif(struct pdb_mol*, short, double*, double*);
/* convert pdb molecular data struct to gen struct */
struct gen_mol *mol_pdb_gen(struct pdb_mol*);
/* convert pdb molecular data struct to xyz struct */
struct xyz_mol *mol_pdb_xyz(struct pdb_mol*);
/* convert pdb struct to zmt molecular data struct */
struct zmt_mol *mol_pdb_zmt(struct pdb_mol*, struct zmt_mol*);

#endif
