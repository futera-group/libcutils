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

#ifndef ZF_LIB_MOL_CUB_H
#define ZF_LIB_MOL_CUB_H

#include <stdio.h>

/* symbolic constants */
#define CUB_CHECK_THRD    0.001  /* threshold for comparison of real numbers */

/* type of incompatibility */
#define CUB_ICMP_NATOM        1  /* number of atoms */
#define CUB_ICMP_ATYPE        2  /* atomic type */
#define CUB_ICMP_GRID         3  /* grid size */
#define CUB_ICMP_LORIG        4  /* origin of lattice */
#define CUB_ICMP_LVEC         5  /* lattice vectors */

/* CUBe atomic data struct */
struct cub_atom {
  unsigned num;             /* atomic number */
  double coord[3];          /* atomic coordinates */
  };

/* CUBe grid paramater struct */
struct cub_grid {
  double origin[3];         /* lattice origin */
  double vector[3][3];      /* lattice vectors */
  unsigned range[3];        /* number of points */
  };

/* CUBE data struct */
struct cub_mol {
  char *title;              /* title line */
  char *desc;               /* description line */
  unsigned n_atoms;         /* number of atom */
  struct cub_atom *atom;    /* atom array */
  struct cub_grid *grid;    /* grid parameters */
  double ***data;           /* data 3D array */
  };

/* allocate new cube molecule struct */
struct cub_mol *mol_cub_new(void);
/* allocate array of cube atom data structures */
struct cub_atom *mol_cub_atom_new(unsigned);
/* allocate new cube grid parameters struct */
struct cub_grid *mol_cub_grid_new(void);
/* allocate new grid data array */
double*** mol_cub_data_new(struct cub_grid*);

/* free memory allocated for cube molecule struct */
void mol_cub_free(struct cub_mol*);
/* free memory allocated for array of atom data structures */
void mol_cub_atom_free(struct cub_atom*);
/* free memory allocated for grid parameter data structure */
void mol_cub_grid_free(struct cub_grid*);
/* free memory allocated for grid data array */
void mol_cub_data_free(struct cub_grid*, double***);

/* create copy of cube grid parameters struct */
struct cub_grid *mol_cub_grid_copy_new(struct cub_grid*);

/* shift data in cube and apply periodic boundary conditions */
void mol_cub_shift(struct cub_mol*, double*);

/* check if two cube data structs are compatible */
short mol_cub_check_compatibility(struct cub_mol*, struct cub_mol*,
  short*, char*, short);

/* read molecular data from cube file */
void mol_cub_read(struct cub_mol*, char*);
/* read molecular data from open cube file */
void mol_cub_fread(struct cub_mol*, FILE*);

/* write molecular data to file in cube format */
void mol_cub_write(struct cub_mol*, char*);
/* write molecular data to file in cube format */
void mol_cub_fwrite(struct cub_mol*, FILE*);
/* write grid parameters to file in cube format */
void mol_cub_fwrite_parm(struct cub_grid*, unsigned, FILE*);

/* convert cub molecular data struct to xyz struct */
struct xyz_mol *mol_cub_xyz(struct cub_mol*);

#endif
