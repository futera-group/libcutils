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

#ifndef ZF_LIB_MOL_ZMT_H
#define ZF_LIB_MOL_ZMT_H

#include <stdio.h>

#define ZMT_FRMT_VAL     0  /* Z-Matrix with values */
#define ZMT_FRMT_VAR     1  /* Z-Matrix with variables */

/* Z-Matrix atom struct */
struct zmt_atom {
  char *name;               /* atomic name (PDB) */
  unsigned num;             /* atomic number */
  unsigned bond_id;         /* bond atom ID */
  double bond_val;          /* bond length */
  unsigned angle_id;        /* angle atom ID */
  double angle_val;         /* angle value */
  unsigned dihed_id;        /* dihedral angle atom ID */
  double dihed_val;         /* dihedral angle */
  };

/* Z-Matrix molecular struct */
struct zmt_mol {
  unsigned n_atoms;         /* number of atoms */
  struct zmt_atom *atom;    /* vector of atoms */
  };

/* Z-Matrix atom struct */
struct zmt_line {
  unsigned num;             /* atomic number */
  unsigned bond_id;         /* bond atom ID */
  char *bond_sym;           /* bond length symbol */
  unsigned angle_id;        /* angle atom ID */
  char *angle_sym;          /* angle symbol */
  unsigned dihed_id;        /* dihedral angle atom ID */
  char *dihed_sym;          /* dihedral angle symbol */
  };

/* Z-Matrix internal-coordinate data structure */
struct zmt_icrd {
  char *sym;                /* internal coordinate symbol */
  double val;               /* internal coordinate value */
  };

/* allocate new zmt molecule struct */
struct zmt_mol *mol_zmt_new(void);
/* allocate new zmt atomic struct */
struct zmt_atom *mol_zmt_atom_new(unsigned);
/* allocate new zmt line data structure */
struct zmt_line *mol_zmt_line_new(void);
/* allocate memory for a new zmt internal-coordinate data structure */
struct zmt_icrd *mol_zmt_icrd_new(void);

/* free memory allocated for zmt molecule struct */
void mol_zmt_free(struct zmt_mol*);
/* free memory allocated for array of zmt atom data structures */
void mol_zmt_atom_free(struct zmt_atom*);
/* free memory allocated for zmt line data structure */
void mol_zmt_line_free(struct zmt_line*);
/* free memory allocated for zmt internal-coordinae data structure */
void mol_zmt_icrd_free(struct zmt_icrd*);

/* copy zmt atom data */
void mol_zmt_atom_copy(struct zmt_atom*, struct zmt_atom*);
/* calculate bond vector for selected atom */
void mol_zmt_atom_vec(double*, double*, double*, double*,
  double, double, double);

/* read molecular data from zmt file */
void mol_zmt_read(struct zmt_mol*, char*);
/* read molecular data from zmt file */
void mol_zmt_fread(struct zmt_mol*, FILE*);

/* write molecular data to file in zmt format */
void mol_zmt_write(struct zmt_mol*, short, char*);
/* write molecular data to file in zmt format */
void mol_zmt_fwrite(struct zmt_mol*, short, FILE*);

/* convert zmt molecular data struct to xyz struct */
struct xyz_mol *mol_zmt_xyz(struct zmt_mol*);

#endif
