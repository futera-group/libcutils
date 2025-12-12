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

#ifndef ZF_LIB_PRG_MOLDEN_H
#define ZF_LIB_PRG_MOLDEN_H

/* symbolic constants */

#define MLD_DATA_ATOM   1  /* atomic coordinates */
#define MLD_DATA_GTO    2  /* Gaussian-type orbitals */
#define MLD_DATA_STO    3  /* Slater-type orbitals */
#define MLD_DATA_MO     4  /* molecular orbitals */
#define MLD_DATA_5D     5  /* spherical d functions */
#define MLD_DATA_5D7F   6  /* spherical d & f functions */
#define MLD_DATA_5D10F  7  /* spherical d and cartesian f functions */
#define MLD_DATA_7F     8  /* spherical f functions */
#define MLD_DATA_9G     9  /* spherical g functions */
#define MLD_DATA_SCF   10  /* SCF convergence */
#define MLD_DATA_GEO   11  /* Geometry-optimization convergence */
#define MLD_DATA_GEOC  12  /* Geometry-optimization structures */
#define MLD_DATA_FREQ  13  /* vibrational frequencies */
#define MLD_DATA_FRC   14  /* vibrational atomic coordinates */
#define MLD_DATA_FRNC  15  /* vibrational normal-mode atomic dispacements */
#define MLD_DATA_FRI   16  /* vibrational intensities */

/* molden atomic data */
struct mld_atom {
  char *sym;                 /* symbol */
  double crd[3];             /* cartesian coordinates */
  unsigned num;              /* atomic numbers */
  };

/* molden data struct - vibrational frequencies */
struct mld_freq {
  double freq;               /* frequency [cm-1] */
  double **displ;            /* coordinate displacement */
  double ir;                 /* IR intensity [KM/mol] */
  };

/* molden data struct */
struct mld_dat {
  unsigned n_atoms;          /* number of atoms */
  unsigned n_vib_modes;      /* number of vibrational modes */
  struct mld_atom *atom;     /* atomic data */
  struct mld_freq *freq;     /* vibrational frequencies */
  };

/* allocate memory for molden data struct */
struct mld_dat *mld_dat_new(void);
/* free memory allocated for molden data struct */
struct mld_dat* mld_dat_free(struct mld_dat*);

/* allocate memory for molden atom data struct */
struct mld_atom *mld_atom_new(unsigned);
/* free memory allocated for molden atom data struct */
struct mld_atom* mld_atom_free(struct mld_atom*, unsigned);

/* allocate array of vibrational-frequency data structs */
struct mld_freq *mld_freq_new(unsigned);
/* free memory allocated for molden frequency data struct array */
struct mld_freq* mld_freq_free(struct mld_freq*, unsigned, unsigned);

/* convert molden data to gaussian file format */
struct gauss_dat* mld_dat_gauss(struct mld_dat*);

/* read data from molden data file */
void mld_read(struct mld_dat*, char*);

/* convert data type symbol to internal ID */
short mld_data_type_id(char*);
/* convert data type internal ID to symbol */
char* mld_data_type_name(short);

/* estimate electric transition dipole moment from the normal mode displ. */
void mld_freq_dip(struct mld_atom*, double**, unsigned, double*);

#endif
